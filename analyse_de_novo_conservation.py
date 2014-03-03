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

from analyse_de_novos import AnalyseDeNovos

class AnalyseDeNovoConservation(AnalyseDeNovos):
    """ class to analyse clustering of de novo events via site specific 
    mutation rates
    """
    
    def analyse_de_novos(self, de_novos, weights):
        """ find the probability of getting de novos with a given mean distance
        
        The probability is the nnumber of simulations where the mean distance
        between simulated de novos is less than the observed distance.
        
        Args:
            de_novos: list of de novos within a gene
            weights: WeightedChoice object to randomly choose positions within
                a gene using site specific mutation rates.
        
        Returns:
            mean distance for the observed de novos and probability of obtaining
            a mean distance less than the observed distance
        """
        
        observed_conservation, sim_prob = "NA", "NA"
        sample_n = len(de_novos)
        if sample_n < 2:
            return (observed_conservation, sim_prob)
        
        dist = self.simulate_distribution(weights, sample_n, self.max_iter)
        
        cds_positions = self.convert_de_novos_to_cds_positions(de_novos)
        observed_conservation = self.get_mean_conservation(cds_positions)
        
        pos = bisect.bisect_right(dist, observed_conservation)
        sim_prob = (1 + pos)/(1 + len(dist))
        
        if type(observed_conservation) != "str":
            observed_conservation = "{0:0.1f}".format(observed_conservation)
        
        return (observed_conservation, sim_prob)
    
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
            
            conservation = self.get_mean_conservation(positions)
            distribution.append(conservation)
        
        distribution = sorted(distribution)
        
        return distribution
    
    def get_mean_conservation(self, positions):
        """ gets the mean distance between two or more CDS positions
        
        Args:
            positions: list of numbers
        
        Returns:
            provides the mean distance of the position pair distances
        """
        
        assert len(positions) > 1
        
        scores = []
        for position in positions:
            score = self.transcript.get_conservation_score(position)
            scores.append(score)
        
        return sum(scores)/len(scores)
    

