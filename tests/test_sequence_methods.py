""" unit test the SequenceMethods class
"""

from __future__ import division
import unittest

from denovonear.transcript import Transcript

class TestSequenceMethodsPy(unittest.TestCase):
    """ unit test the SequenceMethods class
    """
    
    def setUp(self):
        """ construct a Transcript object to add sequence to
        """
        
        self.gene = self.construct_gene()
    
    def construct_gene(self, name='TEST', chrom='1', start=100, end=200,
            strand='+', exons=[(100, 120), (180, 200)],
            cds=[(110, 120), (180, 190)]):
        
        tx = Transcript(name, chrom, start, end, strand)
        tx.set_exons(exons, cds)
        tx.set_cds(cds)
        
        return tx
    
    def test_get_position_on_chrom(self):
        """ test that get_position_on_chrom() works correctly
        """
        
        # check CDS positions for a gene on the forward strand, across two exons
        self.assertEqual(self.gene.get_position_on_chrom(1), 111)
        self.assertEqual(self.gene.get_position_on_chrom(10), 120)
        self.assertEqual(self.gene.get_position_on_chrom(11), 180)
        
        self.assertEqual(self.gene.get_position_on_chrom(10, 1), 121)
        self.assertEqual(self.gene.get_position_on_chrom(11, -1), 179)
        
        # check CDS positions for gene on the reverse strand
        self.gene = self.construct_gene(strand='-')
        self.assertEqual(self.gene.get_position_on_chrom(1), 189)
        self.assertEqual(self.gene.get_position_on_chrom(10), 180)
        self.assertEqual(self.gene.get_position_on_chrom(11), 120)
        
        self.assertEqual(self.gene.get_position_on_chrom(10, -1), 179)
        self.assertEqual(self.gene.get_position_on_chrom(11, 1), 121)
        
        # run through all the CDS positions, and make sure that converting
        # between coordinates gives the same position as original
        for pos in range(self.gene.get_start(), self.gene.get_end()):
            if self.gene.in_coding_region(pos):
                cds = self.gene.chrom_pos_to_cds(pos)
                converted = self.gene.get_position_on_chrom(cds['pos'])
                self.assertEqual(pos, converted)
        
        # now try converting for a gene on the plus strand
        self.gene = self.construct_gene()
        for pos in range(self.gene.get_start(), self.gene.get_end()):
            if self.gene.in_coding_region(pos):
                cds = self.gene.chrom_pos_to_cds(pos)
                converted = self.gene.get_position_on_chrom(cds['pos'])
                self.assertEqual(pos, converted)
    
    def test_get_codon_number_for_cds_position(self):
        """ test that get_codon_number_for_cds_position() works correctly
        """
        
        self.assertEqual(self.gene.get_codon_number_for_cds_position(0), 0)
        self.assertEqual(self.gene.get_codon_number_for_cds_position(1), 0)
        self.assertEqual(self.gene.get_codon_number_for_cds_position(2), 0)
        self.assertEqual(self.gene.get_codon_number_for_cds_position(3), 1)
        self.assertEqual(self.gene.get_codon_number_for_cds_position(4), 1)
        self.assertEqual(self.gene.get_codon_number_for_cds_position(30000), 10000)
    
    def test_get_position_within_codon(self):
        """ test that get_position_within_codon() works correctly
        """
        
        self.assertEqual(self.gene.get_position_within_codon(0), 0)
        self.assertEqual(self.gene.get_position_within_codon(1), 1)
        self.assertEqual(self.gene.get_position_within_codon(2), 2)
        self.assertEqual(self.gene.get_position_within_codon(3), 0)
    
    def test_add_genomic_sequence(self):
        """ test that add_genomic_sequence() works correctly
        """
        
        self.gene = self.construct_gene(start=0, end=10, exons=[(0, 4), (6, 10)],
            cds=[(2, 4), (6, 8)])
        
        gdna = "AGAGGCCTATT"
        # check that we can add gDNA sequence, even if we haven't got any
        # reference CDS to check
        self.gene.add_genomic_sequence(gdna)
        
        self.gene.add_cds_sequence("AGGCTA")
        
        # check that we can add gDNA after adding cDNA
        self.gene.add_genomic_sequence(gdna)
        self.assertEqual(self.gene.get_cds_sequence(), "AGGCTA")
        
        # now check for a gene on the reverse strand
        gdna = "ATAGGCCTCTT"
        self.gene = self.construct_gene(start=0, end=10, exons=[(0, 4), (6, 10)],
            cds=[(2, 4), (6, 8)], strand='-')
        self.gene.add_cds_sequence("AGGCTC")
        self.gene.add_genomic_sequence(gdna)
        self.assertEqual(self.gene.get_cds_sequence(), "AGGCTC")
        
        # check that we raise an error if the CDS predicted from the gDNA
        # doesn't match the CDS already provided
        self.gene.add_cds_sequence("TTTTTT")
        with self.assertRaises(ValueError):
            self.gene.add_genomic_sequence(gdna)
    
    def test_add_genomic_sequence_expanded(self):
        """ test that add_genomic_sequence() works correctly with extra sequence
        """
        
        # now check when the 5' and 3' sequence is extended beyond the gene
        start = 1
        end = 9
        exons = [(1, 4), (6, 9)]
        cds = [(2, 4), (6, 8)]
        self.gene = self.construct_gene(start=start, end=end, exons=exons, cds=cds)
        
        gdna = "AAAGGCCTTT"
        self.gene.add_cds_sequence("AGGCTT")
        
        self.gene.add_genomic_sequence(gdna, offset=1)
        self.assertEqual(self.gene.get_cds_sequence(), "AGGCTT")
    
    def test_add_genomic_sequencE_without_cds_coords(self):
        """ test that error is raised if we add gDNA without CDS coords
        """
        
        a = Transcript("a", "1", 10, 20, "+")
        a.set_exons([(10, 20)], [(10, 20)])
        
        with self.assertRaises(ValueError):
            a.add_genomic_sequence('CGTAGACTGTACGCATCGATT', offset=5)
        
        a.set_cds([(10, 20)])
        a.add_genomic_sequence('CGTAGACTGTACGCATCGATT', offset=5)
    
    def test_get_centered_sequence(self):
        """ test that get_centered_sequence() works correctly
        """
        
        start = 0
        end = 10
        exons = [(0, 4), (6, 10)]
        cds = [(2, 4), (6, 8)]
        self.gene = self.construct_gene(start=start, end=end, exons=exons, cds=cds)
        
        gdna = "AAAGGCCTTT"
        self.gene.add_cds_sequence("AGGCTT")
        self.gene.add_genomic_sequence(gdna, offset=0)
        
        # test CDS positions: start, end, and spanning the exon boundaries
        self.assertEqual(self.gene.get_centered_sequence(2), "AAG")
        self.assertEqual(self.gene.get_centered_sequence(3), "AGG")
        self.assertEqual(self.gene.get_centered_sequence(4), "GGC")
        self.assertEqual(self.gene.get_centered_sequence(6), "CCT")
        self.assertEqual(self.gene.get_centered_sequence(7), "CTT")
        self.assertEqual(self.gene.get_centered_sequence(8), "TTT")
        
        # test that positions outside the CDS raise errors
        with self.assertRaises(ValueError):
            self.gene.get_centered_sequence(-1)
        
        with self.assertRaises(ValueError):
            self.gene.get_centered_sequence(10)
    
    def test_get_centered_sequence_minus_strand(self):
        """ test that get_centered_sequence() works correctly on the minus strand
        """
        
        strand = '-'
        start = 0
        end = 10
        exons = [(0, 4), (6, 10)]
        cds = [(2, 4), (6, 8)]
        self.gene = self.construct_gene(start=start, end=end, exons=exons,
            cds=cds, strand=strand)
        
        gdna = "AAAGGCCTTTT"
        self.gene.add_cds_sequence("AGGCTT")
        self.gene.add_genomic_sequence(gdna, offset=0)
        
        # test CDS positions: start, end, and spanning the exon boundaries
        self.assertEqual(self.gene.get_centered_sequence(2), "AAA")
        self.assertEqual(self.gene.get_centered_sequence(3), "AAG")
    
    def test_get_centered_sequence_long_region(self):
        """ test that get_centered_sequence() works correctly with longer regions
        """
        
        start = 0
        end = 10
        exons = [(0, 4), (6, 10)]
        cds = [(2, 4), (6, 8)]
        self.gene = self.construct_gene(start=start, end=end, exons=exons, cds=cds)
        
        gdna = "AAAGGCCTTT"
        self.gene.add_cds_sequence("AGGCTT")
        self.gene.add_genomic_sequence(gdna, offset=0)
        
        # test CDS positions: start, end, and spanning the exon boundaries
        self.assertEqual(self.gene.get_centered_sequence(2, length=5), "AAAGG")
        
        # check that we raise an error for regions with even numbered lengths
        with self.assertRaises(ValueError):
            self.gene.get_centered_sequence(2, length=4)
    
    def test_get_codon_sequence(self):
        """ test that get_codon_sequence() works correctly
        """
        
        start = 0
        end = 10
        exons = [(0, 4), (6, 10)]
        cds = [(2, 4), (6, 8)]
        self.gene = self.construct_gene(start=start, end=end, exons=exons,
            cds=cds)
        
        gdna = "AAAGGCCTTT"
        self.gene.add_cds_sequence("AGGCTT")
        self.gene.add_genomic_sequence(gdna, offset=0)
        
        # self.gene.start = 0
        # self.gene.end = 10
        # self.gene.exons = [(0, 4), (6, 10)]
        # self.gene.cds = [(2, 4), (6, 8)]
        # self.gene.cds_min = 2
        # self.gene.cds_max = 8
        # self.gene.upstream_sequence = "A"
        # self.gene.cds_sequence = "AGGCTT"
        # self.gene.downstream_sequence = "T"
        
        self.assertEqual(self.gene.get_codon_sequence(0), "AGG")
        self.assertEqual(self.gene.get_codon_sequence(1), "CTT")
        
        # check that codon positions outside the CDS region raise errors
        with self.assertRaises(ValueError):
            self.gene.get_codon_sequence(-1)
        with self.assertRaises(ValueError):
            self.gene.get_codon_sequence(3)
    
    def test_get_codon_sequence_minus_strand(self):
        """ test that get_codon_sequence() works correctly
        """
        
        strand = '-'
        start = 0
        end = 10
        exons = [(0, 4), (6, 10)]
        cds = [(2, 4), (6, 8)]
        self.gene = self.construct_gene(start=start, end=end, exons=exons,
            cds=cds, strand=strand)
        
        gdna = "AAAGGCCTTTT"
        self.gene.add_cds_sequence("AGGCTT")
        self.gene.add_genomic_sequence(gdna, offset=0)
        
        self.assertEqual(self.gene.get_codon_sequence(0), "AGG")
        self.assertEqual(self.gene.get_codon_sequence(1), "CTT")
        
        # check that codon positions outside the CDS region raise errors
        with self.assertRaises(ValueError):
            self.gene.get_codon_sequence(-1)
        with self.assertRaises(ValueError):
            self.gene.get_codon_sequence(3)
    
    def test_translate(self):
        """ test that translate() works correctly
        """
        
        # check a couple of codons, including the stop codon
        self.assertEqual(self.gene.translate("AAA"), "K")
        self.assertEqual(self.gene.translate("AAC"), "N")
        self.assertEqual(self.gene.translate("TAG"), "*")
        
        # raise an error for non-IUPAC base containing codons
        with self.assertRaises(ValueError):
            self.gene.translate("FFF")
        
        # can translate multicodon sequences
        self.assertEqual(self.gene.translate("AAAAAC"), "KN")
        
        # raise an error if we translate sequence with an incomplete codon
        with self.assertRaises(ValueError):
            self.assertEqual(self.gene.translate("AA"), "KN")
            self.assertEqual(self.gene.translate("AAAAA"), "KN")
        
