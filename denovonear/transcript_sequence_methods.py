""" function to translate DNA codons to single character amino acids
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import math

IS_PYTHON2 = sys.version_info[0] == 2

if IS_PYTHON2:
    from string import maketrans
else:
    maketrans = str.maketrans

class SequenceMethods(object):
    
    transdict = maketrans("acgtuACGTU", "tgcaaTGCAA")
    
    def get_position_in_cds(self, chr_position):
        """ figure out the coding position relative to the CDS start site
        """
        
        return self.get_coding_distance(self.get_cds_start(), chr_position)
    
    def get_position_on_chrom(self, cds_position):
        """ figure out the chromosome position of a CDS site
        
        Args:
            cds_position:  position of a variant in CDS
        
        Returns:
            chromosome bp position of the CDS site
        """
        
        # cache the exon boundaries in CDS distances
        if not hasattr(self, "exon_to_cds"):
            self.exon_to_cds = {}
            
            for start, end in self.cds:
                # get the positions of the exon boundaries in CDS distance from
                # the start site
                start_cds = self.get_coding_distance(self.get_cds_start(), start)
                end_cds = self.get_coding_distance(self.get_cds_start(), end)
                
                # cache the CDS positions of the exon boundaries
                self.exon_to_cds[start] = start_cds
                self.exon_to_cds[end] = end_cds
        
        # quickly find the exon containing the CDS position
        for start, end in self.cds:
            start_cds = min(self.exon_to_cds[start], self.exon_to_cds[end])
            end_cds = max(self.exon_to_cds[start], self.exon_to_cds[end])
            
            if start_cds <= cds_position <= end_cds:
                break
        
        # convert the CDS position to a chromosomal coordinate
        if self.get_strand() == "+":
            return start + (cds_position - start_cds)
        else:
            return start + (end_cds - cds_position)
    
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
        
        # if the genomic DNA was for a transcript on the reverse strand, reorient
        # it to the plus strand.
        if self.get_strand() == '-':
            gdna = self.reverse_complement(gdna)
        
        self.gdna_offset = offset
        self.genomic_sequence = gdna
        
        cds_seq = ""
        for start, end in self.cds:
            start_bp = abs(start - self.get_start()) + offset
            end_bp = abs(end - self.get_start()) + 1 + offset
            bases = self.genomic_sequence[start_bp:end_bp]
            if self.strand == "+":
                cds_seq += bases
            else:
                cds_seq = self.reverse_complement(bases) + cds_seq
        
        if cds_seq[0:50] == self.cds_sequence[1:51]:
            self.fix_transcript_off_by_one_bp()
            self.add_genomic_sequence(gdna, offset)
            cds_seq = self.cds_sequence
        
        # do a sanity check to check that we've got the right cds sequence, this
        # fails for at least one gene (CCDC18), which begins with a N, and
        # throws the coordinates off
        if cds_seq != self.cds_sequence:
            raise ValueError("Coding sequence from gene coordinates doesn't"
                "match coding sequence obtained from Ensembl.\nTranscript:"
                "{0}\n{1}\n\nshould be\n{2}\n".format(self.get_name(), cds_seq,
                self.cds_sequence))
    
    def fix_transcript_off_by_one_bp(self):
        """ This fixes ACTL7A, which has the CDS start and end off by a bp.
        """
        
        offset = 1
        if self.strand == "+":
            self.cds[0] = (self.cds[0][0] - offset, self.cds[0][-1])
            self.cds[-1] = (self.cds[-1][0], self.cds[-1][-1] - offset)
            self.exons[0] = (self.exons[0][0] - offset, self.exons[0][-1])
            self.exons[-1] = (self.exons[-1][0], self.exons[-1][-1] - offset)
        else:
            self.cds[0] = (self.cds[0][0] + offset, self.cds[0][-1])
            self.cds[-1] = (self.cds[-1][0], self.cds[-1][-1] + offset)
            self.exons[0] = (self.exons[0][0] + offset, self.exons[0][-1])
            self.exons[-1] = (self.exons[-1][0], self.exons[-1][-1] + offset)
    
    def fix_coding_sequence_length(self):
        """ correct the coding sequence of a transcript, if it misses bases.
        
        Some transcripts don't cover the full length of a gene, they terminate
        in the middle of an amino acid. This is rare, and it is rarer that we
        select these transcripts for analysis. TNS3 is an example, where the
        de novos in the gene can only be contained within a transcript that is
        incomplete at the 3' end. This causes problems when we try to translate
        the final codon.
        
        We simply extend the coding sequence 1-2 bases to make a complete codon.
        """
        
        diff = len(self.cds_sequence) % 3
        end = self.get_cds_end()
        
        if diff != 0:
            diff = 3 - diff
            
            if self.strand == "+":
                self.cds[-1] = (self.cds[-1][0], self.cds[-1][-1] + diff)
                
                start_bp = abs(end - self.get_start()) + self.gdna_offset
                end_bp = abs(end - self.get_start()) + diff + self.gdna_offset
                self.cds_sequence += self.genomic_sequence[start_bp:end_bp]
            elif self.strand == "-":
                self.cds[0] = (self.cds[0][0] - diff, self.cds[0][1])
                start_bp = abs(self.get_end() - end) + self.gdna_offset
                end_bp = abs(self.get_end() - end) + diff + self.gdna_offset
                self.cds_sequence += self.reverse_complement(self.genomic_sequence[start_bp:end_bp])
        
        self.cds_min = int(self.cds[0][0])
        self.cds_max = int(self.cds[-1][-1])
    
    def reverse_complement(self, seq):
        """ reverse complement a DNA or RNA sequence
        """
        
        return str(seq).translate(self.transdict)[::-1]
    
    def get_trinucleotide(self, pos):
        """ obtains the trinucleotide sequence around a position
        
        Args:
            pos: integer for chromosomal nucleotide position.
        
        Returns:
            trinucleotide e.g. 'TCC', centered around the required position,
            given with respect to the plus strand.
        """
        
        assert pos >= 0
        assert pos > self.get_start() - self.gdna_offset and pos < self.get_end() + self.gdna_offset
        
        offset = self.get_start()
        sequence_pos = pos - offset + self.gdna_offset
        tri = self.genomic_sequence[sequence_pos - 1:sequence_pos + 2]
        
        return tri
    
    def get_codon_sequence(self, codon_number):
        """ get the codon sequence for a given codon_number
        """
        
        assert codon_number >= 0
        assert codon_number <= len(self.cds_sequence) / 3
        
        start = codon_number * 3
        end = start + 3
        
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
