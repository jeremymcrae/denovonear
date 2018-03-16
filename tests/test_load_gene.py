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

import unittest
import tempfile
import shutil

from denovonear.load_gene import get_transcript_lengths, construct_gene_object, \
    get_de_novos_in_transcript, get_transcript_ids, load_gene, \
    count_de_novos_per_transcript, minimise_transcripts
from denovonear.transcript import Transcript
from denovonear.ensembl_requester import EnsemblRequest

class TestLoadGenePy(unittest.TestCase):
    """ unit test functions to load genes
    """
    
    @classmethod
    def setUpClass(self):
        self.temp_dir = tempfile.mkdtemp()
        self.ensembl = EnsemblRequest(self.temp_dir, genome_build="grch37")
    
    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.temp_dir)
    
    def set_transcript(self):
        """ construct a transcript for a known gene
        """
        
        exon_ranges=[(120933859, 120934019), (120934219, 120934356),
            (120935876, 120936296)]
        cds_ranges=[(120934225, 120934356), (120935876, 120936013)]
        
        expected = Transcript("ENST00000242577", '12', 120933859, 120936296, "+")
        expected.set_exons(exon_ranges, cds_ranges)
        expected.set_cds(cds_ranges)
        
        cds = "ATGTGCGACCGAAAGGCCGTGATCAAAAATGCGGACATGTCGGAAGAGATGCAACAGGACTC" \
            "GGTGGAGTGCGCTACTCAGGCGCTGGAGAAATACAACATAGAGAAGGACATTGCGGCTCATATC" \
            "AAGAAGGAATTTGACAAGAAGTACAATCCCACCTGGCATTGCATCGTGGGGAGGAACTTCGGTA" \
            "GTTATGTGACACATGAAACCAAACACTTCATCTACTTCTACCTGGGCCAAGTGGCCATTCTTCT" \
            "GTTCAAATCTGGTTAA"
        
        genomic = "GGGCGGGGCCGGCGGGAGCAGGGCGGGGCCTGAGCACTAGGCGGCGGCGGCTGGCGTGGG" \
            "GCTGCTTAGATGCGCCACGGTTTCGGTAGCGACGGTATCTCTAGCCGGGCCTGAGCTGTGCTAGCA" \
            "CCTCCCCCAGGAGACCGTTGCAGTCGGCCAGCCCCCTTCTCCACGGTGAGAAACTCGGGGGGCCAG" \
            "GGGGTGTCCTCGCTGCCTTATTTCGCCCCACTCCGGACTTAGCCCTCCGCGTAGCCCGCGCTTCCT" \
            "GAGAAGTGGGGTGGGGGGCGTCGTCCCGTGGTGGCGCCGGCCGGGGTGGGGGCAGTTAGTGCCTGG" \
            "GGGGCGCGGCCCAACTCAACCCCTTACCCCAGGCCTTGCCCACTAGGTAACCATGTGCGACCGAAA" \
            "GGCCGTGATCAAAAATGCGGACATGTCGGAAGAGATGCAACAGGACTCGGTGGAGTGCGCTACTCA" \
            "GGCGCTGGAGAAATACAACATAGAGAAGGACATTGCGGCTCATATCAAGAAGGTGAGGATGGGCGC" \
            "GGGGGCCGATACGCAGCCGGGAGCAGGGGGTTCCTTCCCCCCGATCCTGCTTTCCTAAGGGCGCCT" \
            "GACAGGTCCCGGGAATACTGCTGGCGGCTTGGGGCGTAGAAGCTTCCAGAAAGGACGCAGATGCAT" \
            "TTTGCGCTCCTGTGGAGAAGACCAGACCCCCGGCGTCCGAAGTTTTTTTTTTTTTTTTTTTAATTA" \
            "CCCAGCTCCGCGGGGGGAAAGCGCCACCTAGCAACGGTATCTAAGATCAGGGAGCAGCGGTTCCCC" \
            "CTTCTGTGTGGTTCCTGCGCCGAGGATCCATCTGGGTGTTCCGGAGGGGGGAGCTGCGTGGGTGTT" \
            "TCCAGCCGGGCCGGGAGGAGATCTTGCCAGCCTTCCAGTGGGGAGTTGAGGGAAGGTGGTGGGTGG" \
            "TGGCGGGGCTGGGGGCTGGGGTAGGGGCTTGGTAAATGGCAGTCTAGAAAGCCGGCAGGACTGCCA" \
            "ACTTCTCGAGCAGTGTTTGCTGGAAGGGAAGAAAGCTGGCAGCCTAAGCCGTGGGAGGGTTCCAGT" \
            "CGAGAATGGGAAGATGAAAGACTTCAGATGGAACAGAAATAAATGCCTTTTTTGACAAACGCAGCA" \
            "GTGCGTGCCTCTAGCTTGCAAGAGCGTTACTCCCCTTCATAGCTTTAAAAGGTTTTCGCACTGCGT" \
            "GCAGTTAGAGTAGCTAAATCTTGTGTGACGCTCCACAAACACTTGTAAGAATTTTGCAGAGAAAGA" \
            "TAACCGTTGCCACCCAATGCCCCCCACAGGCATTCTACTCCCCAGTACCTCTTAGGGTGGGAGAAA" \
            "TGGTGAAGAGTTGTTCCTACAACTTGCTAACCTAGTGGACAGGGTAGTAGATTAGCATCATCCGGA" \
            "TAGATGTGAAGAGGACGGCTGTTTGGATAATAATTAAGGATAAAATTTGGCCAGTTGACAGATTCT" \
            "GTTTCCAGCAGTTTTTACAGCAACAGTGGAGTGCTTCAGTATTGTGTTCCTGTAAATTTAATTTTG" \
            "ATCCGCAATCATTTGGTATACAATGCTGTTTGAAGTTTTGTCCTATTGGAAAAGTCTTGTGTTGCA" \
            "GGGGTGCAGTTAAGATCTTTGTGATGAGGAATGGGATGGGCTAATTTTTTGCCGTTTTCTTGGAAT" \
            "TGGGGGCATGGCAAATACAGTAGGGTAGTTTAGTTCTCTACACAGAACATGATAAACTACACCTGT" \
            "TGATGTCACCGTCTGTCAATGAATATTATAGAAGGTATGAAGGTGTAATTACCATAATAACAAAAC" \
            "ACCCTGTCTTTAGGGCTGACCTTTCGTCCTTTGACCTCCTCAGCCTCCATTCCCATCTTCGCTCAG" \
            "ACTGCAAGTATGTTTGTATTAATGTACTATGTAGGCGGCTTGGAGCTGGGGAACATTCTTTCATTC" \
            "TAAGAATTTGCAGATGCTGACGTTCCTCCTTTCTGCCCCTACAGGCTCTGGCTTATCCAAGAGGCA" \
            "AACACTGACCTCTGGTAATTAAAATCCTAGTTCTTTTCTTTTGTCTTTTCCAGGAATTTGACAAGA" \
            "AGTACAATCCCACCTGGCATTGCATCGTGGGGAGGAACTTCGGTAGTTATGTGACACATGAAACCA" \
            "AACACTTCATCTACTTCTACCTGGGCCAAGTGGCCATTCTTCTGTTCAAATCTGGTTAAAAGCATG" \
            "GACTGTGCCACACACCCAGTGATCCATCCAAAAACAAGGACTGCAGCCTAAATTCCAAATACCAGA" \
            "GACTGAAATTTTCAGCCTTGCTAAGGGAACATCTCGATGTTTGAACCTTTGTTGTGTTTTGTACAG" \
            "GGCATTCTCTGTACTAGTTTGTCGTGGTTATAAAACAATTAGCAGAATAGCCTACATTTGTATTTA" \
            "TTTTCTATTCCATACTTCTGCCCACGTTGTTTTCTCTCAAAATCCATTCCTTTAAAAAATAAATCT" \
            "GATGCAGATGTGTATGTGTGTG"
        
        expected.add_cds_sequence(cds)
        expected.add_genomic_sequence(genomic, offset=10)
        
        return expected
    
    def test_get_transcript_lengths(self):
        """ check that get_transcript_lengths() works correctly
        """
        
        transcript_ids = ["ENST00000392509", "ENST00000242577", "ENST00000549649"]
        lengths = {'ENST00000242577': 89, 'ENST00000392509': 89, 'ENST00000549649': 42}
        
        self.assertEqual(get_transcript_lengths(self.ensembl, transcript_ids), lengths)
    
    def test_construct_gene_object(self):
        """
        """
        
        transcript_id = "ENST00000242577"
        transcript = construct_gene_object(self.ensembl, transcript_id)
        
        expected = self.set_transcript()
        
        self.assertEqual(transcript, expected)
        self.assertEqual(transcript.get_genomic_sequence(), expected.get_genomic_sequence())
        self.assertEqual(transcript.get_cds_sequence(), expected.get_cds_sequence())
    
    def test_get_de_novos_in_transcript(self):
        """ test that we can identify de novos within the CDS of a transcript
        """
        
        exon_ranges = [(10, 20), (30, 40), (90, 100)]
        cds_ranges = [(30, 40), (90, 95)]
        
        # define a simple transcript
        tx = Transcript("test1", '1', 10, 100, "+")
        tx.set_exons(exon_ranges, cds_ranges)
        tx.set_cds(cds_ranges)
        
        # check that only the site in the CDS is returned
        sites = [15, 35, 100]
        self.assertEqual(get_de_novos_in_transcript(tx, sites), [35])
        
        # check that we can return multiple sites in the CDS
        sites = [15, 35, 90]
        self.assertEqual(get_de_novos_in_transcript(tx, sites), [35, 90])
        
        # check if we pass in an empty list, we get one back
        self.assertEqual(get_de_novos_in_transcript(tx, []), [])
    
    def test_get_transcript_ids(self):
        """
        """
        
        hgnc = "DYNLL1"
        lengths = get_transcript_ids(self.ensembl, hgnc)
        
        expected = {'ENST00000548342': 89,
            'ENST00000549989': 89, 'ENST00000392509': 89,
            'ENST00000392508': 89, 'ENST00000242577': 89,
            'ENST00000550845': 67, 'ENST00000550178': 67,
            'ENST00000548214': 67, 'ENST00000552870': 47,
            'ENST00000549649': 42}
        
        self.assertEqual(lengths, expected)
    
    def test_load_gene(self):
        """ check that we correctly load the suitable transcripts for a gene
        """
        
        # define a hgnc symbol to load transcript for, and de novo sites to
        # check against
        hgnc = "DYNLL1"
        sites = [120934226, 120936012]
        transcripts = load_gene(self.ensembl, hgnc, sites)
        
        # define the expected transcript, and make sure that it is in the list
        # of suitable transcripts. There can be multiple transcripts return if
        # more than one transcript of the maximal length includes all de novos.
        expected = self.set_transcript()
        self.assertIn(expected, transcripts)
        
        # and make sure if none of the de novos fall in a suitable transcript,
        # then we raise an error.
        sites = [100, 200]
        with self.assertRaises(IndexError):
            transcripts = load_gene(self.ensembl, hgnc, sites)
    
    def test_count_de_novos_per_transcript(self):
        """ test that we count de novos in transcripts correctly
        """
        
        hgnc = "DYNLL1"
        sites = [120934226, 120936012]
        counts = count_de_novos_per_transcript(self.ensembl, hgnc, sites)
        
        expected = {'ENST00000549649': {'len': 42, 'n': 1},
            'ENST00000548214': {'len': 67, 'n': 1},
            'ENST00000242577': {'len': 89, 'n': 2},
            'ENST00000392508': {'len': 89, 'n': 2},
            'ENST00000392509': {'len': 89, 'n': 2},
            'ENST00000549989': {'len': 89, 'n': 2},
            'ENST00000550178': {'len': 67, 'n': 1},
            'ENST00000552870': {'len': 47, 'n': 1},
            'ENST00000550845': {'len': 67, 'n': 1},
            'ENST00000548342': {'len': 89, 'n': 2}}
        
        self.assertEqual(counts, expected)
        
        # TODO: add test case for error from gene where no protein coding
        # TODO: transcript is available
        
    def test_minimise_transcripts(self):
        """ test that minimise_transcripts() works correctly
        """
        
        # run through a test case for a single gene
        hgnc = "DYNLL1"
        sites = [120934226, 120936012]
        counts = minimise_transcripts(self.ensembl, hgnc, sites)
        expected = {'ENST00000242577': {'len': 89, 'n': 2},
            'ENST00000392508': {'len': 89, 'n': 2},
            'ENST00000392509': {'len': 89, 'n': 2},
            'ENST00000549989': {'len': 89, 'n': 2},
            'ENST00000548342': {'len': 89, 'n': 2}}
        
        self.assertEqual(counts, expected)
        
        # check that when we don't have any de novos, we return an empty list
        self.assertEqual(minimise_transcripts(self.ensembl, hgnc, []), {})
        
        # check that when none of the de novos are in a transcript, we return
        # an empty list.
        self.assertEqual(minimise_transcripts(self.ensembl, hgnc, [100]), {})
