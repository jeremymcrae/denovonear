""" class to analyse clustering of known de novos in genes according to their
distances apart within the gene, and compare that to simulated de novo events
within the same gene.
"""

from __future__ import division
from __future__ import print_function

import itertools

from denovonear.analyse_de_novos import AnalyseDeNovos

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
    
    def get_geometric_mean(self, pos_pairs):
        """ get the geometric mean distance between pair positions
        """
        
        distances = []
        for pos_1, pos_2 in pos_pairs:
            distance = abs(pos_1 - pos_2)
            distances.append(distance)
        
        return self.geomean(distances)
