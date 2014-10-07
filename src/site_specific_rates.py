""" get weight gene mutation rates
"""

from __future__ import print_function

from src.weighted_choice import WeightedChoice


class SiteRates(object):
    """ class to build weighted choice random samplers for nonsense, missense, 
    and functional classes of variants, using site specific mutation rates
    
    Only include the site specific probability if it mutates to a different
    amino acid, or occurs close to an intron/exon boundary, using consequences
    defined at: http://www.ensembl.org/info/genome/variation/predicted_data.html
    """
    
    def __init__(self, gene, mut_dict, masked_sites=None):
        
        self.gene = gene
        self.mut_dict = mut_dict
        self.masked = masked_sites
    
    def get_missense_rates_for_gene(self):
        return self.build_weighted_site_rates_for_gene(self.missense_check)
        
    def get_nonsense_rates_for_gene(self):
        return self.build_weighted_site_rates_for_gene(self.nonsense_check)
    
    def get_functional_rates_for_gene(self):
        return self.build_weighted_site_rates_for_gene(self.functional_check)
    
    def get_synonymous_rates_for_gene(self):
        return self.build_weighted_site_rates_for_gene(self.synonymous_check)
        
    def get_lof_rates_for_gene(self):
        return self.build_weighted_site_rates_for_gene(self.loss_of_function_check)
    
    def get_splice_loss_of_function_rates_for_gene(self):
        return self.build_weighted_site_rates_for_gene(self.splice_lof_check)
    
    def get_splice_region_rates_for_gene(self):
        return self.build_weighted_site_rates_for_gene(self.splice_region_check)
    
    def get_missense_and_splice_region_rates_for_gene(self):
        return self.build_weighted_site_rates_for_gene(self.missense_and_splice_region_check)
    
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
        
        return self.loss_of_function_check(initial_aa, mutated_aa, position) or \
            self.missense_and_splice_region_check(initial_aa, mutated_aa, position)
    
    def splice_lof_check(self, initial_aa, mutated_aa, position):
        """ checks if a variant has a splice_donor or splice_acceptor consequence
        
        These variants are defined as being the two intronic base-pairs adjacent
        to the intron/exon boundary.
        """
        
        return (not self.gene.in_coding_region(position)) and self.boundary_dist < 3
    
    def nonsense_check(self, initial_aa, mutated_aa, position):
        """ checks if two amino acids are a nonsense (eg stop_gained) mutation
        """
        
        return initial_aa != "*" and mutated_aa == "*"
    
    def loss_of_function_check(self, initial_aa, mutated_aa, position):
        """ checks if a site gives a loss of function mutation (ie nonsense or
        affects splice consensus)
        """
        
        return self.nonsense_check(initial_aa, mutated_aa, position) or \
            self.splice_lof_check(initial_aa, mutated_aa, position)
    
    def missense_and_splice_region_check(self, initial_aa, mutated_aa, position):
        """ check for missense, or splice_region_variant (but not lof)
        """
        
        # trim out nonsense mutations such as stop_gained mutations, and splice 
        # site mutations
        if self.loss_of_function_check(initial_aa, mutated_aa, position):
            return False
        
        return self.missense_check(initial_aa, mutated_aa, position) or \
            self.splice_region_check(initial_aa, mutated_aa, position)
    
    def missense_check(self, initial_aa, mutated_aa, position):
        """ checks if two amino acids are a missense mutation (but not nonsense)
        """
        
        # trim out nonsense mutations such as stop_gained mutations, and splice 
        # site mutations
        if self.loss_of_function_check(initial_aa, mutated_aa, position):
            return False
        
        # include the site if it mutates to a different amino acid.
        return initial_aa != mutated_aa 
    
    def splice_region_check(self, initial_aa, mutated_aa, position):
        """ checks if a variant has a splice_donor or splice_acceptor consequence
        """
        
        # catch splice region variants within the exon, and in the appropriate 
        # region of the intron (note that loss of function splice_donor and 
        # splice_acceptor variants have been excluded when we trimmed nonsense).
        if self.gene.in_coding_region(position):
            # check for splice_region_variant inside exon
            return self.boundary_dist < 4
        else:
            # check for splice_region_variant inside intron
            return self.boundary_dist < 9
    
    def synonymous_check(self, initial_aa, mutated_aa, position):
        """ checks if two amino acids are synonymous
        """
        
        return not self.functional_check(initial_aa, mutated_aa, position)
    
    def get_gene_range(self):
        """ get the lowest and highest base positions of a genes CDS
        """
        
        boundary_1 = self.gene.get_cds_start()
        boundary_2 = self.gene.get_cds_end()
        
        start = min(boundary_1, boundary_2)
        end = max(boundary_1, boundary_2)
        
        if self.gene.strand == "+":
            start -= 1
        
        end += 1
        
        return (start, end)
    
    def build_weighted_site_rates_for_gene(self, mut_type):
        """ build a list of sites in a gene that can give missense mutations,  
        along with their weighted probability of the mutation occuring.
        
        Args:
            mut_type: function to check if a site matches the mutation type
        
        Returns:
            WeightedChoice object, for random sampling by the weight of the 
            probabilities
        """
        
        probs = []
        range_start, range_end = self.get_gene_range()
        
        for bp in range(range_start, range_end):
            # ignore sites within masked regions (typically masked because the
            # site has been picked up on alternative transcript)
            if self.masked is not None and self.masked.in_coding_region(bp):
                continue
            
            # ignore sites outside the CDS region
            if bp < min(self.gene.get_cds_start(), self.gene.get_cds_end()) or \
                bp > max(self.gene.get_cds_start(), self.gene.get_cds_end()):
                continue
            
            # get the distances to the closest exon boundaries
            exon_start, exon_end = self.gene.find_closest_exon(bp)
            self.boundary_dist = min(abs(exon_start - bp), abs(exon_end - bp))
            
            # ignore positions outside the exons that are too distant from an
            # intron/exon boundary
            if not self.gene.in_coding_region(bp) and self.boundary_dist >= 9:
                continue
            
            # use a default NA amino acid, which will be changed if the current
            # position lies within the coding sequence
            initial_aa = "NA"
            self.seq = self.gene.get_trinucleotide(bp)
            if self.gene.in_coding_region(bp):
                cds_pos = self.gene.get_position_in_cds(bp)
                codon = self.gene.get_codon_number_for_cds_position(cds_pos)
                self.codon_pos = self.gene.get_position_within_codon(cds_pos)
                self.codon = self.gene.get_codon_sequence(codon)
                initial_aa = self.gene.translate_codon(self.codon)
            else:
                # for non-exonic positions, use the nearest exon boundary
                site = [x for x in [exon_start, exon_end] if abs(x - bp) == self.boundary_dist]
                cds_pos = self.gene.get_position_in_cds(site[0])
            
            # drop the initial base, since we want to mutate to other bases
            bases = ["A", "C", "G", "T"]
            bases.remove(self.seq[1])
            
            for base in bases:
                mutated_aa = initial_aa
                mutated_seq = self.seq[0] + base + self.seq[2]
                if self.gene.in_coding_region(bp):
                    mutated_aa = self.get_mutated_aa(base, self.codon, self.codon_pos)
                
                if mut_type(initial_aa, mutated_aa, bp):
                    probs.append([cds_pos, self.mut_dict[self.seq][mutated_seq]])
        
        return WeightedChoice(probs)


