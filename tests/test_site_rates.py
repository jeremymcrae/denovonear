"""
Copyright (c) 2015 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from __future__ import print_function

import os
import unittest
import itertools

from denovonear.site_specific_rates import get_gene_range, get_mutated_aa, SiteRates
from denovonear.weights import WeightedChoice
from denovonear.transcript import Transcript

class TestSiteRatesPy(unittest.TestCase):
    """ unit test the SiteRates class
    """
    
    bases = set(["A", "T", "G", "C"])
    
    rates = []
    for initial in itertools.product("".join(bases), repeat=3):
        trinuc = "".join(initial)
        
        for base in bases - set(initial[1]):
            changed = "{}{}{}".format(initial[0], base, initial[2])
            rate = '5e-7'
            
            rates.append([trinuc, changed,  '5e-7'])
    
    def setUp(self):
        """ construct a Transcript object to add sequence to
        """
        
        self.transcript = self.construct_gene()
        
        cds = "ATGTCCATAACCAAAGCCTGA"
        genomic = "CCTCCAGATTCACGGGAAGCATGTCCATAAGTAGGGAGATATTTGGTGCTCTCATTTG" \
            "TGGAGACTCTAGCCAAAGCCTGAGTCATGCGTACCATAGATAG"
        
        self.transcript.add_cds_sequence(cds)
        self.transcript.add_genomic_sequence(genomic, offset=10)
        
        self.weights = SiteRates(self.transcript, self.rates)
    
    def construct_gene(self, name='TEST', chrom='1', start=100, end=179,
            strand='+', exons=[(100, 119), (160, 179)],
            cds=[(110, 119), (160, 170)]):
        
        tx = Transcript(name, chrom, start, end, strand)
        tx.set_exons(exons, cds)
        tx.set_cds(cds)
        
        return tx
    
    def test_get_boundary_distance(self):
        """ check the function to get distances to the nearest intron/exon boundary
        """
        
        # check a site upstream of the gene
        self.assertEqual(self.transcript.get_boundary_distance(50), 50)
        
        # check a site at the start of a gene
        self.assertEqual(self.transcript.get_boundary_distance(100), 0)
        
        # check some sites within the first exon
        self.assertEqual(self.transcript.get_boundary_distance(110), 10)
        self.assertEqual(self.transcript.get_boundary_distance(115), 5)
        
        # check sites in the first intron
        self.assertEqual(self.transcript.get_boundary_distance(125), 6)
        self.assertEqual(self.transcript.get_boundary_distance(140), 20)
        
        # check a site in the first exon, as it becomes closer to the next intron
        self.assertEqual(self.transcript.get_boundary_distance(141), 19)
        
        # check a site downstream of the gene
        self.assertEqual(self.transcript.get_boundary_distance(200), 21)
    
    def test_get_codon_info(self):
        """ check the function that checks the codon information for a position
        """
        
        # make sure a site well outside the gene raises an error
        with self.assertRaises(ValueError):
            self.transcript.get_codon_info(50)
        
        # a position near the start site, but upstream of the CDS will raise a
        # different error
        with self.assertRaises(RuntimeError):
            self.transcript.get_codon_info(95)
        
        # check the first base of the CDS
        self.assertEqual(self.transcript.get_codon_info(110),
            {'cds_pos': 0, 'codon_seq': 'ATG', 'intra_codon': 0,
                "codon_number": 0, 'initial_aa': 'M'})
        
        # check the second base of the CDS
        self.assertEqual(self.transcript.get_codon_info(111),
            {'cds_pos': 1, 'codon_seq': 'ATG', 'intra_codon': 1,
                "codon_number": 0, 'initial_aa': 'M'})
        
        # check the third base of the CDS
        self.assertEqual(self.transcript.get_codon_info(112),
            {'cds_pos': 2, 'codon_seq': 'ATG', 'intra_codon': 2,
                "codon_number": 0, 'initial_aa': 'M'})
        
        # check the fourth base of the CDS
        self.assertEqual(self.transcript.get_codon_info(113),
            {'cds_pos': 3, 'codon_seq': 'TCC', 'intra_codon': 0,
                "codon_number": 1, 'initial_aa': 'S'})
        
        # check a site 2 bp into the first intron. We assign this as the
        # position of the closest exon boundary, but without any codon info
        self.assertEqual(self.transcript.get_codon_info(122),
            {'cds_pos': 9, 'codon_seq': None, 'intra_codon': None,
                "codon_number": None, 'initial_aa': None})
    
    def test_site_rates_weights(self):
        """ check the cumulative mutation rates for each consequence type.
        """
        
        wts = self.weights
        
        # these cumulative mutation rates should be correct for each consequence
        # type
        self.assertAlmostEqual(wts["missense"].get_summed_rate(), 2.45e-05, places=7)
        self.assertAlmostEqual(wts["nonsense"].get_summed_rate(), 5e-07, places=7)
        self.assertAlmostEqual(wts["synonymous"].get_summed_rate(), 4e-06, places=7)
        self.assertAlmostEqual(wts["loss_of_function"].get_summed_rate(), 6.5e-06, places=7)
        self.assertAlmostEqual(wts["splice_lof"].get_summed_rate(), 6e-06, places=7)
        self.assertAlmostEqual(wts["splice_region"].get_summed_rate(), 2.05e-05, places=7)
    
    def test_site_rates_sampled(self):
        """ check the sites sampled for each consequence group.
        
        Repeatedly sample sites within the transcript. Given enough samples, we
        should saturate the transcript, and so we can check that we have sampled
        all of those sites.
        """
        
        wts = self.weights
        n = 10000
        
        # there are numerous sites for mutating to a missense change
        self.assertEqual(set([ wts["missense"].choice() for x in range(n) ]),
            set([0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 18, 19, 20]))
        
        # the transcript only has one position where we can get a stop_gained
        self.assertEqual(set([ wts["nonsense"].choice() for x in range(n) ]),
            set([12]))
        
        # the transcript only has a few positions where we can get synonymous changes
        self.assertEqual(set([ wts["synonymous"].choice() for x in range(n) ]),
            set([5, 14, 17, 19]))
        
        # the transcript only has one splice donor site and one splice acceptor
        # site. The splice lof sites are shifted to the nearest exon coordinate
        # in order to be able to check CDS proximity.
        self.assertEqual(set([ wts["splice_lof"].choice() for x in range(n) ]),
            set([9, 10]))
        
        # loss-of-function sites are the union of nonsense and splice_lof sites
        self.assertEqual(set([ wts["loss_of_function"].choice() for x in range(n) ]),
            set([9, 10, 12]))
        
        # splice region variants can occur at the two bp inside the exon, or
        # within 8 bp of the intron/exon boundary (but beyond the splice lof
        # positions). The intronic positions are shifted to the nearest exon
        # coordinate for CDS proximity checking.
        self.assertEqual(set([ wts["splice_region"].choice() for x in range(n) ]),
            set([8, 9, 10, 11]))
    
    def test_get_mutated_aa(self):
        """ check that mutating a codon gives the expected amino acids
        """
        
        # check some codon mutations, including a stop mutation
        self.assertEqual(get_mutated_aa(self.transcript, "C", "AAA", 2), "N")
        self.assertEqual(get_mutated_aa(self.transcript, "A", "TGG", 2), "*")
        
        # a codon mutated to itself gives the expected amino acid
        self.assertEqual(get_mutated_aa(self.transcript, "A", "AAA", 2), "K")
        
        # non-DNA codons raise errors
        with self.assertRaises(ValueError):
            get_mutated_aa(self.transcript, "C", "RRR", 2)
    
    def test_splice_lof_check(self):
        """ check that splice_lof_check() works correctly
        """
        
        # check a site within the intron, but only 2 bp away from the
        # intron/exon boundary
        self.weights.check_position(121)
        self.assertTrue(self.weights.splice_lof_check('', '', 121))
        
        # check a site within the exon, and also 2 bp away from the intron/exon
        # boundary
        self.weights.check_position(117)
        self.assertFalse(self.weights.splice_lof_check('', '', 117))
        
        # check a intron site outside the splice lof positions
        self.weights.check_position(122)
        self.assertFalse(self.weights.splice_lof_check('', '', 122))
    
    def test_nonsense_check(self):
        """ check that nonsense_check() works correctly
        """
        
        self.weights.check_position(161)
        self.assertTrue(self.weights.nonsense_check("N", "*", 161))
        self.assertFalse(self.weights.nonsense_check("N", "G", 161))
        self.assertFalse(self.weights.nonsense_check("*", "G", 161))
        self.assertFalse(self.weights.nonsense_check("*", "*", 161))
    
    def test_missense_check(self):
        """ check that missense_check() works correctly
        """
        
        # missense mutations can either be where the amino acids differ, or a
        # stop site changes to coding an amino acid (technically these are
        # stop_lost, but they carry a missense-like severity).
        self.weights.check_position(161)
        self.assertTrue(self.weights.missense_check("N", "G", 161))
        self.assertTrue(self.weights.missense_check("*", "G", 161))
        
        # don't include stop gained mutations, or stop to stop
        self.assertFalse(self.weights.missense_check("N", "*", 161))
        self.assertFalse(self.weights.missense_check("*", "*", 161))
        
        # the case below shouldn't occur, where the site is in the intron and
        # within the splice lof positions, but somehow the initial and modified
        # amino acids differ, but check this anyway.
        # self.weights.boundary_dist = 2
        self.weights.check_position(121)
        self.assertFalse(self.weights.missense_check("N", "G", 121))
    
    def test_splice_region_check(self):
        """ check that splice_region_check() works correctly
        """
        
        self.weights.check_position(123)
        self.assertTrue(self.weights.splice_region_check("", "", 123))
        
        # check a site in the splice lof positions
        self.weights.check_position(121)
        self.assertFalse(self.weights.splice_region_check("", "", 121))
        
        # check a site just inside the splice region positions
        self.weights.check_position(127)
        self.assertTrue(self.weights.splice_region_check("", "", 128))
        
        # check a site just outside the splice region positions
        self.weights.check_position(130)
        self.assertFalse(self.weights.splice_region_check("", "", 130))
        
        # check an exonic splice region position
        self.weights.check_position(116)
        self.assertTrue(self.weights.splice_region_check("N", "N", 116))
        
        # check an exonic site just beyond the splice region positions
        self.weights.check_position(115)
        self.assertFalse(self.weights.splice_region_check("N", "N", 115))
        
        # check an exonic site inside the splice region positions, but with
        # mutated amino acids aren't splice region variants
        self.weights.check_position(116)
        self.assertFalse(self.weights.splice_region_check("N", "K", 116))
    
    def test_synonymous_check(self):
        """ check that synonymous_check() works correctly
        """
        
        self.weights.check_position(115)
        self.assertTrue(self.weights.synonymous_check("N", "N", 115))
        self.assertTrue(self.weights.synonymous_check("*", "*", 115))
        
        # amino acid changes aren't synonymous
        self.assertFalse(self.weights.synonymous_check("N", "*", 115))
        self.assertFalse(self.weights.synonymous_check("*", "N", 115))
        
        # sites in splice region or splice lof aren't synonymous
        self.weights.check_position(117)
        self.assertFalse(self.weights.synonymous_check("N", "N", 117))
        self.weights.check_position(121)
        self.assertFalse(self.weights.synonymous_check("", "", 121))
    
    def test_get_gene_range(self):
        """ check that get_gene_range() works correctly
        """
        
        self.assertEqual(get_gene_range(self.transcript), {'start': 110, 'end': 170})
    
    def test_check_position_missense_only(self):
        """ check that check_position() works correctly for missense changes
        """
        
        self.weights.clear()
        self.weights.check_position(110)
        
        # mutating the first base of the start codon all has three missense
        # possibilities, and if all those have a uniform rate of 5e-7, then the
        # missense rate should sum to 1.5e-7
        self.assertEqual(self.weights["missense"].get_summed_rate(), 1.5e-6)
        self.assertEqual(self.weights["nonsense"].get_summed_rate(), 0)
        self.assertEqual(self.weights["synonymous"].get_summed_rate(), 0)
        self.assertEqual(self.weights["splice_lof"].get_summed_rate(), 0)
        self.assertEqual(self.weights["splice_region"].get_summed_rate(), 0)
        self.assertEqual(self.weights["loss_of_function"].get_summed_rate(), 0)
        
    def test_check_position_mixed_nonsense(self):
        """ check that check_position() works correctly for mixed changes
        """
        
        self.weights.clear()
        self.weights.check_position(162)
        
        # check a position where one alternate base gives a nonsense change
        self.assertEqual(self.weights["missense"].get_summed_rate(), 1.0e-6)
        self.assertEqual(self.weights["nonsense"].get_summed_rate(), 0.5e-6)
        self.assertEqual(self.weights["synonymous"].get_summed_rate(), 0)
        self.assertEqual(self.weights["splice_lof"].get_summed_rate(), 0)
        self.assertEqual(self.weights["splice_region"].get_summed_rate(), 0)
        self.assertEqual(self.weights["loss_of_function"].get_summed_rate(), 0.5e-6)
    
    def test_check_position_noncoding(self):
        """ check that check_position() works correctly for noncoding changes
        """
        
        self.weights.clear()
        categories = ['missense', 'nonsense', 'synonymous', 'splice_lof', 'splice_region']
        
        # a deep intronic site won't alter the summed rates
        self.weights.check_position(140)
        for x in categories:
            self.assertEqual(self.weights[x].get_summed_rate(), 0)
        
        # an upstream site won't alter the summed rates
        self.weights.check_position(self.transcript.get_start() - 1)
        for x in categories:
            self.assertEqual(self.weights[x].get_summed_rate(), 0)
        
        # a downstream site won't alter the summed rates
        self.weights.check_position(self.transcript.get_end() + 1)
        for x in categories:
            self.assertEqual(self.weights[x].get_summed_rate(), 0)
        
        # a splice lof only affects certain rates
        self.weights.check_position(121)
        self.assertEqual(self.weights["missense"].get_summed_rate(), 0)
        self.assertEqual(self.weights["nonsense"].get_summed_rate(), 0)
        self.assertEqual(self.weights["synonymous"].get_summed_rate(), 0)
        self.assertEqual(self.weights["splice_lof"].get_summed_rate(), 1.5e-6)
        self.assertEqual(self.weights["splice_region"].get_summed_rate(), 0)
        self.assertEqual(self.weights["loss_of_function"].get_summed_rate(), 1.5e-6)
        
        # a splice region only affects certain rates
        self.weights.clear()
        self.weights.check_position(125)
        self.assertEqual(self.weights["missense"].get_summed_rate(), 0)
        self.assertEqual(self.weights["nonsense"].get_summed_rate(), 0)
        self.assertEqual(self.weights["synonymous"].get_summed_rate(), 0)
        self.assertEqual(self.weights["splice_lof"].get_summed_rate(), 0)
        self.assertEqual(self.weights["splice_region"].get_summed_rate(), 1.5e-6)
        self.assertEqual(self.weights["loss_of_function"].get_summed_rate(), 0)
