""" class to analyse clustering of known de novos in genes according to their 
distances apart within the gene, and compare that to simulated de novo events 
within the same gene.
"""

from __future__ import division
from __future__ import print_function

import bisect
import itertools
import math
import operator


class AnalyseDeNovos(object):
    """ class to analyse clustering of de novo events via site specific 
    mutation rates
    """
    def __init__(self, transcript, site_weights, iterations):
        """ initialise the class
        """
        
        self.transcript = transcript
        self.site_weights = site_weights
        self.max_iter = iterations
    
    def analyse_missense(self, de_novo_events):
        """ analyse clustering of missense de novos
        """
        
        weights = self.site_weights.get_missense_rates_for_gene()
        return self.analyse_de_novos(de_novo_events, weights)
    
    def analyse_nonsense(self, de_novo_events):
        """ analyse clustering of nonsense de novos
        """
        
        weights = self.site_weights.get_nonsense_rates_for_gene()
        return self.analyse_de_novos(de_novo_events, weights)
    
    def analyse_functional(self, de_novo_events):
        """ analyse clustering of functional (missense and nonsense) de novos
        """
        
        weights = self.site_weights.get_functional_rates_for_gene()
        return self.analyse_de_novos(de_novo_events, weights)
    
    def analyse_de_novos(self, de_novos, weights):
        """ find the probability of getting de novos with a mean conservation
        
        The probability is the number of simulations where the mean conservation
        between simulated de novos is less than the observed conservation.
        
        Args:
            de_novos: list of de novos within a gene
            weights: WeightedChoice object to randomly choose positions within
                a gene using site specific mutation rates.
        
        Returns:
            mean conservation for the observed de novos and probability of 
            obtaining a mean conservation less than the observed conservation
        """
        
        observed_value, sim_prob = "NA", "NA"
        sample_n = len(de_novos)
        if sample_n < 2:
            return (observed_value, sim_prob)
        
        dist = self.simulate_distribution(weights, sample_n, self.max_iter)
        
        cds_positions = self.convert_de_novos_to_cds_positions(de_novos)
        observed_value = self.get_score(cds_positions)
        
        pos = bisect.bisect_right(dist, observed_value)
        sim_prob = (1 + pos)/(1 + len(dist))
        
        if type(observed_value) != "str":
            observed_value = "{0:0.1f}".format(observed_value)
        
        return (observed_value, sim_prob)
    
    def simulate_distribution(self, weights, sample_n=2, max_iter=100):
        """ creates a distribution of mutation scores in a single gene
        
        Args:
            weights: WeightedChoice object
            sample_n: number of de novo mutations to sample
            max_iter: number of iterations/simulations to run
        """
        
        distribution = []
        iteration = 0
        while iteration < max_iter:
            iteration += 1
            
            positions = []
            while len(positions) < sample_n:
                site = weights.choice()
                positions.append(site)
            
            # the following line is class specific - can do distance clustering,
            # conservation scores
            value = self.get_score(positions)
            distribution.append(value)
        
        distribution = sorted(distribution)
        
        return distribution
    
    def convert_de_novos_to_cds_positions(self, de_novos):
        """ convert cds positions for de novo events into cds positions
        
        Args:
            de_novos: list of chrom bp positions within the transcript
        
        Returns:
            list of positions converted to CDS positions within the transcript
        """
        
        # need to convert the de novo event positions into CDS positions
        cds_start = self.transcript.get_cds_start()
        if self.transcript.strand == "-":
            cds_start = self.transcript.get_cds_end()
        cds_positions = []
        for pos in de_novos:
            pos += 1 # offset to zero based chrom
            try:
                dist = self.transcript.get_coding_distance(cds_start, pos)
            except AssertionError:
                # catch the splice site functional mutations
                (start, end) = self.transcript.find_closest_exon(pos)
                
                start_dist = abs(start - pos)
                end_dist = abs(end - pos)
                
                # if the var is outside the exon, but affects a splice site, 
                # swap it to using the splice site location
                if start_dist < 3:
                    dist = self.transcript.get_coding_distance(cds_start, start)
                elif end_dist < 3:
                    dist = self.transcript.get_coding_distance(cds_start, end)
                else:
                    raise ValueError("distance to exon (" + str(max(start_dist,\
                        end_dist)) + ") > 2 bp for " + str(pos) + " in " + \
                        "transcript " + self.transcript.get_name())
                
            cds_positions.append(dist)
        
        return cds_positions
    
    def geomean(self, values):
        """ get the geometric mean of a list of values
        """
        
        # get the geometric mean, but be careful around values of 0, since
        # without correction, the mean distance would be zero
        if 0 in values:
            # allow for 0s in a geometric mean by shifting everything up one, 
            # then dropping the mean by one at the end
            values = [x + 1 for x in values]
            total = self.product(values)
            mean = total ** (1/len(values))
            mean -= 1
        else:
            total = self.product(values)
            mean = total ** (1/len(values))
        
        return mean


