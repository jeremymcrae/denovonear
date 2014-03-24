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
    """ class to analyse conservation of de novo events via site specific 
    mutation rates
    """
    
    def get_score(self, positions):
        """ gets the mean conservation of two or more CDS positions
        
        Args:
            positions: list of numbers
        
        Returns:
            provides the mean conservation at the gene positions
        """
        
        assert len(positions) > 1
        
        scores = []
        for position in positions:
            score = self.transcript.get_conservation_score(position)
            scores.append(score)
        
        # get the mean score,  but make the score negative, since the simulation
        # looks for a minimal value
        return -(sum(scores)/len(scores))
    

