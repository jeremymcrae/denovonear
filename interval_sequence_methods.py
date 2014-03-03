""" function to translate DNA codons to single character amino acids
"""

from __future__ import division
from __future__ import print_function

import string
import math

class SequenceMethods(object):
    
    def get_position_in_cds(self, chr_position):
        """ figure out the coding position relative to the CDS start site
        """
        
        if self.strand == "+":
            return self.get_coding_distance(self.get_cds_start(), chr_position)
        elif self.strand == "-":
            return self.get_coding_distance(self.get_cds_end(), chr_position)
        else:
            raise ValueError("unknown strand type" + self.strand)
    
    def get_position_on_chrom(self, cds_position):
        """ figure out the chromosome position of a CDS site
        """
        
        for start, end in self.cds:
            if self.strand == "+":
                start_cds = self.get_coding_distance(self.get_cds_start(), start)
                end_cds = self.get_coding_distance(self.get_cds_start(), end)
            else:
                start_cds = self.get_coding_distance(self.get_cds_end(), end)
                end_cds = self.get_coding_distance(self.get_cds_end(), start)
            
            # print(cds_position, start_cds, end_cds, start, end)
            
            if start_cds <= cds_position <= end_cds:
                break
        
        if self.strand == "+":
            return start + (cds_position - start_cds)
        else:
            return end - (cds_position - start_cds)
    
    def get_codon_number_for_cds_position(self, cds_position):
        """ figure out the codon position for a position
        """
        
        return int(math.floor(cds_position / 3))
    
    def get_position_within_codon(self, cds_position):
        """ get the position within a codon (in 0 based format eg 0, 1, 2)
        """
        
        return (cds_position) % 3
    
    def add_cds_sequence(self, cds_dna):
        """ add the CDS sequence
        """
        
        self.cds_sequence = cds_dna
    
    def add_genomic_sequence(self, gdna, offset=0):
        """ add and process the genomic sequence into a CDS sequence
        """
        
        gdna = gdna[offset:]
        self.genomic_sequence = gdna
        
        cds_seq = ""
        for start, end in self.cds:
            if self.strand == "+":
                start_bp = abs(start - self.get_start())
                end_bp = abs(end - self.get_start())
                cds_seq += gdna[start_bp:end_bp + 1]
            else:
                start_bp = abs(self.get_end() - end)
                end_bp = abs(self.get_end() - start)
                cds_seq = gdna[start_bp - 1:end_bp] + cds_seq
        
        # do a sanity check to check that we've got the right cds sequence
        if cds_seq != self.cds_sequence:
            raise ValueError("haven't obtained the right CDS for " + \
                self.get_name() + "\n" + cds_seq)
        
        self.cds_sequence = cds_seq
        
        if self.strand == "+":
            up_pos = abs(self.get_cds_start() - self.get_start())
            down_pos = abs(self.get_cds_end() - self.get_start())
        else:
            up_pos = abs(self.get_cds_end() - self.get_end())
            down_pos = abs(self.get_cds_start() - self.get_end())
        
        self.upstream_sequence = gdna[up_pos]
        self.downstream_sequence = gdna[down_pos]
        
    def reverse_complement(self, seq):
        """ reverse complement a DNA or RNA sequence
        """
        
        transdict = string.maketrans("acgtuACGTU", "tgcaaTGCAA")
        
        return seq.translate(transdict)[::-1]
    
    def get_trinucleotide_around_cds_position(self, cds_position):
        """ obtains the trinucleotide sequence around a cds position
        """
        
        if cds_position == 0:
            tri = self.upstream_sequence + self.cds_sequence[:2]
        elif cds_position == len(self.cds_sequence) - 1:
            tri = self.cds_sequence[-2:] + self.downstream_sequence
        else:
            tri = self.cds_sequence[cds_position - 1:cds_position + 2]
        
        return tri
    
    def get_codon_sequence(self, codon_number):
        """ get the codon sequence for a given codon_number
        """
        
        start = codon_number * 3
        end = start + 3
        # end = codon_number * 3
        # start = end - 3
        codon = self.cds_sequence[start:end]
        
        return codon
    
    def translate_codon(self, codon_sequence):
        """ translate a DNA codon to a single character amino acid
        """
        
        t = {"AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N", \
            "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T", \
            "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S", \
            "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I", \
            "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H", \
            "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P", \
            "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R", \
            "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L", \
            "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D", \
            "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A", \
            "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G", \
            "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V", \
            "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "Y", \
            "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S", \
            "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C", \
            "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"}
        
        codon_sequence = codon_sequence.upper()
        
        return t[codon_sequence]
    
        
        
    
