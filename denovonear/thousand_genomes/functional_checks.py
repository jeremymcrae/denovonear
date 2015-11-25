""" script to obtain 1000 genomes variation data in CDS exons
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import


class FunctionalChecks(object):
    """ intended to be a subclass of Extract1000Genomes, where this contains
    functions to check if a variant has any functional consequences
    """
    
    def get_functional_status(self, var):
        """ figure out if a variant is lof, missense (functional), or synonymous
        """
        
        pos = var.pos
        
        # get the distances to the closest exon boundaries
        exon_start, exon_end = self.gene.find_closest_exon(pos)
        pos_dist = min(abs(exon_start - pos), abs(exon_end - pos))
        ref_end = pos + len(var.ref) - 1
        alt_end = pos + len(var.alt[0]) - 1
        ref_end_dist = abs(exon_start - ref_end), abs(exon_end - ref_end)
        alt_end_dist = abs(exon_start - alt_end), abs(exon_end - alt_end)
        self.boundary_dist = min(pos_dist, ref_end_dist, alt_end_dist)
        
        (ref, alt) = self._get_major_and_minor_allele_by_freq(var)
        
        (lof_check, missense_check) = self.check_exon_spanning_deletion(pos, ref, exon_start, exon_end)
        if lof_check:
            return "lof"
        elif missense_check:
            return "missense"
        elif self.check_lof_indel(ref, alt):
            return "lof"
        elif self.check_missense_indel(ref, alt):
            return "missense"
        
        (lof_check, missense_check) = self.check_snv(pos, ref, alt)
        if lof_check:
            return "lof"
        elif missense_check:
            return "missense"
        else:
            return "synonymous"
    
    def check_lof_indel(self, ref, alt):
        """ check if a variants is a frameshift indel
        """
        
        lof = False
        if len(alt) > 1 and (len(alt) - 1) % 3 != 0:
            lof = True
        elif len(ref) > 1 and (len(ref) - 1) % 3 != 0:
            lof = True
        
        return lof
    
    def check_missense_indel(self, ref, alt):
        """ check if a variant is an in-frame indel
        """
        
        missense = False
        if len(alt) > 1 and (len(alt) - 1) % 3 == 0:
            missense = True
        elif len(ref) > 1 and (len(ref) - 1) % 3 == 0:
            missense = True
        
        return missense
    
    def check_exon_spanning_deletion(self, pos, ref, exon_start, exon_end):
        """ check for deletions that span the exon boundaries
        
        I won't consider insertions, since insertions upstream of the exon
        probably wouldn't impact the transcript, whereas insertions that
        extend in distance past the exon end would simply extend the
        transcript sequence, so it only matters if the insertion is in
        frame or not. Obviously this will not always be true, eg
        insertions that introduce new intron/exon boundaries.
        """
        
        lof = False
        missense = False
        
        # don't check SNV variants in this function, only deletions
        if len(ref) == 1:
            return (lof, missense)
        
        if pos < exon_start or pos + len(ref) > exon_end:
            ref_end = pos + len(ref) - 1
            # check for deletions that surround an intron/exon boundary
            if pos < exon_start and ref_end > exon_start:
                lof = True
            elif pos < exon_end and ref_end > exon_end:
                lof = True
            # check for deletions that finish just short of an exon
            elif pos < exon_start and ref_end < exon_start:
                if abs(ref_end - exon_start) < 3:
                    lof = True
                elif abs(ref_end - exon_start) < 9:
                    missense = True
            # check for deletions that start just after an exon
            elif pos > exon_end:
                if self.boundary_dist < 3:
                    lof = False
                elif self.boundary_dist < 9:
                    missense = True
        
        return (lof, missense)
    
    def check_snv(self, pos, ref, alt):
        """ get the functional status of a SNV
        """
        
        if self.gene.strand == "-":
            ref = self.gene.reverse_complement(ref)
            alt = self.gene.reverse_complement(alt)
            pos -= 1
        
        # get the codon containing the variant, and the intra-codon position
        cds_pos = self.gene.convert_chr_pos_to_cds_positions(pos)
        codon_number = self.gene.get_codon_number_for_cds_position(cds_pos)
        codon_pos = self.gene.get_position_within_codon(cds_pos)
        codon = self.gene.get_codon_sequence(codon_number)
        
        # make sure the reference codon contains the reference allele (the
        # alleles might have been swapped since I select the more frequent
        # allele as the reference, rather than directly from the reference
        # sequence)
        if codon[codon_pos] == alt:
            codon = list(codon)
            codon[codon_pos] = ref
            codon = "".join(codon)
        
        # figure out if the variant impacts the amino acid.
        # Occasionally the variant occurs in an incomplete transcript (eg
        # rs28394186 in ENST00000336769). Those may give odd length codons
        # which we shall ignore (since they are so rare).
        initial_aa = self.gene.translate_codon(codon)
        mutated_aa = initial_aa
        if self.gene.in_coding_region(pos):
            mutated_aa = self.get_mutated_aa(alt, codon, codon_pos)
        
        pos += 1
        exon_start, exon_end = self.gene.find_closest_exon(pos)
        self.boundary_dist = min(abs(exon_start - pos), abs(exon_end - pos))
        if self.gene.in_coding_region(pos):
            self.boundary_dist += 1
        
        lof = False
        missense = False
        if self.missense_check(initial_aa, mutated_aa, pos) or \
                self.splice_region_check(initial_aa, mutated_aa, pos):
            missense = True
        elif self.nonsense_check(initial_aa, mutated_aa, pos) or \
                self.splice_lof_check(initial_aa, mutated_aa, pos):
            lof = True
        
        return (lof, missense)
    
    
    
    
    
