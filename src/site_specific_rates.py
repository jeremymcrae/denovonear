""" get weight gene mutation rates
"""

from __future__ import print_function

from src.weighted_choice import WeightedChoice


class SiteRates(object):
    """ class to build weighted choice random samplers for nonsense, missense, 
    and functional classes of variants, using site specific mutation rates
    """
    
    def __init__(self, gene, mut_dict):
        
        self.gene = gene
        self.mut_dict = mut_dict
    
    def get_missense_rates_for_gene(self):
        return self.build_weighted_site_rates_for_gene(self.missense_check)
    
    def get_nonsense_rates_for_gene(self):
        return self.build_weighted_site_rates_for_gene(self.nonsense_check)
    
    def get_functional_rates_for_gene(self):
        return self.build_weighted_site_rates_for_gene(self.functional_check)
    
    def get_mutated_aa(self, base, codon, codon_pos):
        """ find the amino acid resulting from a base change to a codon
        """
          
        # figure out what the mutated codon is
        mutated_codon = list(codon)
        mutated_codon[codon_pos] = base
        mutated_codon = "".join(mutated_codon)
        
        return self.gene.translate_codon(mutated_codon)
    
    def functional_check(self, initial_aa, mutated_aa, position):
        """ checks if two amino acids are have a functional mutation (either
        nonsense or missense)
        """
        
        # only include the site specific probability if it mutates to a 
        # different amino acid, or occurs close to an intron/exon boundary
        # (using the splice_region_variant definition at: 
        # http://www.ensembl.org/info/genome/variation/predicted_data.html)
        if initial_aa != mutated_aa or self.exon_start_dist < 4 or \
                self.exon_end_dist < 4:
            return True
        
        return False
    
    def nonsense_check(self, initial_aa, mutated_aa, position):
        """ checks if two amino acids are a nonsense (eg stop_gained) mutation
        """
        
        # print("stop", initial_aa != "*" and mutated_aa == "*")
        # print("near splice", not self.gene.in_coding_region(position) and \
        #     (self.exon_start_dist < 3 or self.exon_end_dist < 3))
        # print("not in coding", not self.gene.in_coding_region(position))
        # print("distance less than 3 bp", self.exon_start_dist < 3 or self.exon_end_dist < 3)
        
        return initial_aa != "*" and mutated_aa == "*" or \
            (not self.gene.in_coding_region(position) and \
            (self.exon_start_dist < 3 or self.exon_end_dist < 3))
    
    def missense_check(self, initial_aa, mutated_aa, position):
        """ checks if two amino acids are a missense mutation (but not nonsense)
        """
        
        # trim out nonsense mutations such as stop_gained mutations, and splice 
        # site mutations
        if initial_aa != "*" and mutated_aa == "*" or \
            (not self.gene.in_coding_region(position) and \
            (self.exon_start_dist < 3 or self.exon_end_dist < 3)):
            return False
        
        # if initial_aa != "*" and mutated_aa == "*" or \
        #     self.exon_start_dist < 3 and self.exon_end_dist < 3:
        #     return False
        
        # only include the site specific probability if it mutates to a 
        # different amino acid, or occurs close to an intron/exon boundary.
        if initial_aa != mutated_aa:
            return True
        
        return False   
    
    def get_gene_range(self):
        """ get the start and end points for a gene
        """
        
        start = self.gene.get_cds_start()
        if self.gene.strand == "+":
            start -= 1
        
        end = self.gene.get_cds_end() + 1
        
        return (start, end)
    
    def build_weighted_site_rates_for_gene(self, mut_type):
        """ build a list of sites in a gene that can give missense mutations, along 
        with their weighted probability of the mutation occuring.
        
        Args:
            mut_type: function to check if a site matches the mutation type
        
        Returns:
            WeightedChoice object, for random sampling by the weight of the mutation
            probabilities
        """
        
        probs = []
        range_start, range_end = self.get_gene_range()
        
        for bp in range(range_start, range_end):
            if not self.gene.in_coding_region(bp):
                continue
            
            cds_pos = self.gene.get_position_in_cds(bp)
            
            codon_number = self.gene.get_codon_number_for_cds_position(cds_pos)
            self.codon_pos = self.gene.get_position_within_codon(cds_pos)
            self.codon = self.gene.get_codon_sequence(codon_number)
            
            # get the distances to the closest exon boundaries
            exon_start, exon_end = self.gene.find_closest_exon(bp)
            self.exon_start_dist = abs(exon_start - bp)
            self.exon_end_dist = abs(exon_end - bp)
            
            self.seq = self.gene.get_trinucleotide_around_cds_position(cds_pos)
            
            # drop the initial base, since we want to mutate to other bases
            bases = ["A", "C", "G", "T"]
            bases.remove(self.seq[1])
            
            initial_aa = self.gene.translate_codon(self.codon)
            for base in sorted(bases):
                mutated_aa = self.get_mutated_aa(base, self.codon, self.codon_pos)
                mutated_seq = self.seq[0] + base + self.seq[2]
                if mut_type(initial_aa, mutated_aa, bp):
                    probs.append([cds_pos, self.mut_dict[self.seq][mutated_seq]])
        
        return WeightedChoice(probs)


