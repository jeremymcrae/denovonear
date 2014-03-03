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

class AnalyseDeNovoClustering(AnalyseDeNovos):
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
        
        observed_distance, sim_prob = "NA", "NA"
        sample_n = len(de_novos)
        if sample_n < 2:
            return (observed_distance, sim_prob)
        
        dist = self.simulate_distribution(weights, sample_n, self.max_iter)
        
        cds_positions = self.convert_de_novos_to_cds_positions(de_novos)
        observed_distance = self.get_mean_distance_between_positions(cds_positions)
        
        pos = bisect.bisect_right(dist, observed_distance)
        sim_prob = (1 + pos)/(1 + len(dist))
        
        if type(observed_distance) != "str":
            observed_distance = "{0:0.1f}".format(observed_distance)
        
        return (observed_distance, sim_prob)
    
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
            
            distance = self.get_mean_distance_between_positions(positions)
            distribution.append(distance)
        
        distribution = sorted(distribution)
        
        return distribution
    
    def get_mean_distance_between_positions(self, positions):
        """ gets the mean distance between two or more CDS positions
        
        Args:
            positions: list of numbers
        
        Returns:
            provides the mean distance of the position pair distances
        """
        
        assert len(positions) > 1
        
        if len(positions) == 2:
            return abs(positions[0] - positions[1])
        
        pos_pairs = itertools.combinations(positions, 2)
        
        return self.get_geometric_mean(pos_pairs)
        # return self.get_arithmetic_mean(pos_pairs)
        # return self.get_mean_from_closest_neighbors(positions)
    
    def get_geometric_mean(self, pos_pairs):
        """ get the geometric mean distance between pair positions
        """
        
        distances = []
        for pos_1, pos_2 in pos_pairs:
            distance = abs(pos_1 - pos_2)
            distances.append(distance)
        
        return self.geomean(distances)
        
    def get_arithmetic_mean(self, pos_pairs):
        """ get the arithmetic mean distance between pair positions
        """
        
        distances = []
        for pos_1, pos_2 in pos_pairs:
            distance = abs(pos_1 - pos_2)
            distances.append(distance)
            
        return sum(distances)/len(distances)
    
    def get_mean_from_closest_neighbors(self, positions):
        """ find the mean distance to the closest neighbor for each variant
        """
        
        run_pairs = set([])
        
        distances = []
        for pos_1 in range(len(positions)):
            bp_1 = positions[pos_1]
            
            min_distance = 999999999
            pair = ("NA", "NA")
            for pos_2 in range(len(positions)):
                # ignore if we are using the same list position
                if pos_1 == pos_2:
                    continue
                
                bp_2 = positions[pos_2]
                distance = abs(bp_1 - bp_2)
                if distance < min_distance:
                    min_distance = distance
                    pair = tuple(sorted([pos_1, pos_2]))
            
            # only include the distance if the pair has not been included  
            # previously (this avoids duplicates of variants that match to each 
            # other)
            if pair not in run_pairs:
                distances.append(distance)
                run_pairs.add(pair)
        
        return self.geomean(distances)
    
    def product(self, iterable):
        return reduce(operator.mul, iterable, 1)

