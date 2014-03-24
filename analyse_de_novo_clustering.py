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
    
    def get_score(self, positions):
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
    

