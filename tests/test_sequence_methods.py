""" unit test the SequenceMethods class
"""

from __future__ import division
import unittest

from denovonear.transcript_sequence import SequenceMethods
from denovonear.transcript import Transcript

class TestSequenceMethodsPy(unittest.TestCase):
    """ unit test the SequenceMethods class
    """
    
    def setUp(self):
        """ construct a Transcript object to add sequence to
        """
        
        chrom = "1"
        start = 100
        end = 200
        name = "TEST"
        strand = "+"
        exons = [(100, 120), (180, 200)]
        cds = [(110, 120), (180, 190)]
        
        self.gene = Transcript(name, start, end, strand, chrom, exons, cds)
    
    def test_get_position_on_chrom(self):
        """ test that get_position_on_chrom() works correctly
        """
        
        # check CDS positions for a gene on the forward strand, across two exons
        self.assertEqual(self.gene.get_position_on_chrom(1), 111)
        self.assertEqual(self.gene.get_position_on_chrom(10), 120)
        self.assertEqual(self.gene.get_position_on_chrom(11), 180)
        
        # check CDS positions for gene on the reverse strand
        self.gene.strand = "-"
        del self.gene.exon_to_cds
        self.assertEqual(self.gene.get_position_on_chrom(1), 189)
        self.assertEqual(self.gene.get_position_on_chrom(10), 180)
        self.assertEqual(self.gene.get_position_on_chrom(11), 120)
        
        # run through all the CDS positions, and make sure that converting
        # between coordinates gives the same position as original
        for pos in range(self.gene.get_start(), self.gene.get_end()):
            if self.gene.in_coding_region(pos):
                cds = self.gene.chrom_pos_to_cds(pos)
                converted = self.gene.get_position_on_chrom(cds)
                self.assertEqual(pos, converted)
        
        # now try converting for a gene on the plus strand
        self.gene.strand = '+'
        del self.gene.exon_to_cds
        for pos in range(self.gene.get_start(), self.gene.get_end()):
            if self.gene.in_coding_region(pos):
                cds = self.gene.chrom_pos_to_cds(pos)
                converted = self.gene.get_position_on_chrom(cds)
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
        
        self.gene.start = 0
        self.gene.end = 10
        self.gene.exons = [(0, 4), (6, 10)]
        self.gene.cds = [(2, 4), (6, 8)]
        self.gene.cds_min = 2
        self.gene.cds_max = 8
        
        # check that we get an error if we haven't got any reference CDS to check
        with self.assertRaises(AttributeError):
            self.gene.add_genomic_sequence("AAA")
        
        gdna = "AGAGGCCTAT"
        self.gene.cds_sequence = "AGGCTA"
        
        self.gene.add_genomic_sequence(gdna)
        self.assertEqual(self.gene.cds_sequence, "AGGCTA")
        
        # now check for a gene on the reverse strand
        self.gene.strand = "-"
        self.gene.cds_sequence = "GAGCCT"
        self.gene.add_genomic_sequence(gdna)
        self.assertEqual(self.gene.cds_sequence, "GAGCCT")
        
        self.gene.cds_sequence = "TTTTTT"
        with self.assertRaises(ValueError):
            self.gene.add_genomic_sequence(gdna)
    
    def test_add_genomic_sequence_off_by_one(self):
        """ test that add_genomic_sequence() works when the seq is off by one
        """
        
        self.gene.start = 0
        self.gene.end = 70
        self.gene.exons = [(5, 60)]
        self.gene.cds = [(5, 60)]
        self.gene.cds_min = 6
        self.gene.cds_max = 61
        
        self.gene.cds_sequence =           "ATGTGGGCTCCACCAGCAGCAATCATGGGGGATGGGCCCACCAAGAAGGTGGGCAA"
        self.gene.add_genomic_sequence("GGGGATGTGGGCTCCACCAGCAGCAATCATGGGGGATGGGCCCACCAAGAAGGTGGGCAACCAGGCCCC")
        self.assertEqual(self.gene.cds_sequence, "ATGTGGGCTCCACCAGCAGCAATCATGGGGGATGGGCCCACCAAGAAGGTGGGCAA")
    
    def test_fix_transcript_off_by_one_bp(self):
        """ test that fix_transcript_off_by_one_bp is correct
        """
        
        self.gene.start = 0
        self.gene.end = 70
        self.gene.exons = [(5, 60)]
        self.gene.cds = [(5, 60)]
        self.gene.cds_min = 6
        self.gene.cds_max = 61
        
        self.gene.fix_transcript_off_by_one_bp()
        self.assertEqual(self.gene.cds, [(4, 59)])
        
        self.gene.strand = "-"
        self.gene.fix_transcript_off_by_one_bp()
        self.assertEqual(self.gene.cds, [(5, 60)])
        
    
    def test_add_genomic_sequence_expanded(self):
        """ test that add_genomic_sequence() works correctly with extra sequence
        """
        
        # now check when the 5' and 3' sequence is extended beyond the gene
        self.gene.start = 1
        self.gene.end = 9
        self.gene.exons = [(1, 4), (6, 9)]
        self.gene.cds = [(2, 4), (6, 8)]
        self.gene.cds_min = 2
        self.gene.cds_max = 8
        
        gdna = "AAAGGCCTTT"
        self.gene.cds_sequence = "AGGCTT"
        
        self.gene.add_genomic_sequence(gdna, offset=1)
        self.assertEqual(self.gene.cds_sequence, "AGGCTT")
    
    def test_get_trinucleotide(self):
        """ test that get_trinucleotide() works correctly
        """
        
        self.gene.start = 0
        self.gene.end = 10
        self.gene.exons = [(0, 4), (6, 10)]
        self.gene.cds = [(2, 4), (6, 8)]
        self.gene.cds_min = 2
        self.gene.cds_max = 8
        self.gene.gdna_offset = 0
        self.gene.genomic_sequence = "AAAGGCCTTT"
        
        # test CDS positions: start, end, and spanning the exon boundaries
        self.assertEqual(self.gene.get_trinucleotide(2), "AAG")
        self.assertEqual(self.gene.get_trinucleotide(3), "AGG")
        self.assertEqual(self.gene.get_trinucleotide(4), "GGC")
        self.assertEqual(self.gene.get_trinucleotide(6), "CCT")
        self.assertEqual(self.gene.get_trinucleotide(7), "CTT")
        self.assertEqual(self.gene.get_trinucleotide(8), "TTT")
        
        # test that positions outside the CDS raise errors
        with self.assertRaises(AssertionError):
            self.gene.get_trinucleotide(-1)
            
        with self.assertRaises(AssertionError):
            self.gene.get_trinucleotide(10)
        
        # test when we define the sequence by the add_genomic_sequence method
        self.gene.start = 1
        self.gene.end = 9
        self.gene.exons = [(1, 4), (6, 9)]
        self.gene.cds = [(2, 4), (6, 8)]
        self.gene.cds_sequence = "AGGCTT"
        self.gene.add_genomic_sequence("AAAGGCCTTT", offset=1)
        self.assertEqual(self.gene.get_trinucleotide(2), "AAG")
    
    def test_get_trinucleotide_minus_strand(self):
        """ test that get_trinucleotide() works correctly on the minus strand
        """
        
        self.gene.strand = "-"
        self.gene.start = 0
        self.gene.end = 10
        self.gene.exons = [(0, 4), (6, 10)]
        self.gene.cds = [(2, 4), (6, 8)]
        self.gene.cds_min = 2
        self.gene.cds_max = 8
        self.gene.gdna_offset = 0
        self.gene.genomic_sequence = "AAAGGCCTTT"
        
        # test CDS positions: start, end, and spanning the exon boundaries
        self.assertEqual(self.gene.get_trinucleotide(2), "AAG")
        self.assertEqual(self.gene.get_trinucleotide(3), "AGG")
    
    def test_get_codon_sequence(self):
        """ test that get_codon_sequence() works correctly
        """
        
        self.gene.start = 0
        self.gene.end = 10
        self.gene.exons = [(0, 4), (6, 10)]
        self.gene.cds = [(2, 4), (6, 8)]
        self.gene.cds_min = 2
        self.gene.cds_max = 8
        self.gene.upstream_sequence = "A"
        self.gene.cds_sequence = "AGGCTT"
        self.gene.downstream_sequence = "T"
        
        self.assertEqual(self.gene.get_codon_sequence(0), "AGG")
        self.assertEqual(self.gene.get_codon_sequence(1), "CTT")
        
        # check that codon positions outside the CDS region raise errors
        with self.assertRaises(AssertionError):
            self.gene.get_codon_sequence(-1)
        with self.assertRaises(AssertionError):
            self.gene.get_codon_sequence(3)
    
    def test_get_codon_sequence_minus_strand(self):
        """ test that get_codon_sequence() works correctly
        """
        
        self.gene.strand = '-'
        self.gene.start = 0
        self.gene.end = 10
        self.gene.exons = [(0, 4), (6, 10)]
        self.gene.cds = [(2, 4), (6, 8)]
        self.gene.cds_min = 2
        self.gene.cds_max = 8
        self.gene.upstream_sequence = "A"
        self.gene.cds_sequence = "AGGCTT"
        self.gene.downstream_sequence = "T"
        
        self.assertEqual(self.gene.get_codon_sequence(0), "AGG")
        self.assertEqual(self.gene.get_codon_sequence(1), "CTT")
        
        # check that codon positions outside the CDS region raise errors
        with self.assertRaises(AssertionError):
            self.gene.get_codon_sequence(-1)
        with self.assertRaises(AssertionError):
            self.gene.get_codon_sequence(3)
    
    def test_translate(self):
        """ test that translate() works correctly
        """
        
        # check a couple of codons, including the stop codon
        self.assertEqual(self.gene.translate("AAA"), "K")
        self.assertEqual(self.gene.translate("AAC"), "N")
        self.assertEqual(self.gene.translate("TAG"), "*")
        
        # raise an error for non-IUPAC base containing codons
        with self.assertRaises(KeyError):
            self.gene.translate("FFF")
        
        # can translate multicodon sequences
        self.assertEqual(self.gene.translate("AAAAAC"), "KN")
        
        # raise an error if we translate sequence with an incomplete codon
        with self.assertRaises(KeyError):
            self.assertEqual(self.gene.translate("AA"), "KN")
            self.assertEqual(self.gene.translate("AAAAA"), "KN")
        
