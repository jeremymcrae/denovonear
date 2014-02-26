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
        
        dist = self.build_distance_distribution(weights, sample_n, self.max_iter)
        
        cds_positions = self.convert_de_novos_to_cds_positions(de_novos)
        observed_distance = self.get_mean_distance_between_positions(cds_positions)
        
        pos = bisect.bisect_right(dist, observed_distance)
        sim_prob = pos/len(dist)
        
        if type(observed_distance) != "str":
            observed_distance = "{0:0.1f}".format(observed_distance)
        
        return (observed_distance, sim_prob)
    
    def build_distance_distribution(self, weights, sample_n=2, max_iter=100):
        """ creates a distribution of distances between mutations in a single gene
        
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
    
    def convert_de_novos_to_cds_positions(self, de_novos):
        """ convert cds positions for de novo events into cds positions
        
        Args:
            de_novos: list of chrom bp positions within the transcript
        
        Returns:
            list of positions converted to CDS positions within the transcript
        """
        
        # need to convert the de novo event positions into CDS positions
        cds_start = self.transcript.get_cds_start()
        cds_positions = []
        for pos in de_novos:
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
        
        distances = [1 if x==0 else x for x in distances]
        total = self.product(distances)
        mean = total ** (1/len(distances))
        
        return mean
        
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
        
        # get the geometric mean, but be careful around distances of 0, since
        # without correction, the mean distance would be zero
        if 0 in distances:
            # allow for 0s in a geometric mean by shifting everything up one, 
            # then dropping the mean by one at the end
            distances = [x + 1 for x in distances]
            total = self.product(distances)
            mean = total ** (1/len(distances))
            mean -= 1
        else:
            total = self.product(distances)
            mean = total ** (1/len(distances))
        
        return mean
    
    def product(self, iterable):
        return reduce(operator.mul, iterable, 1)



