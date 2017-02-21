"""
Copyright (c) 2015 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import math
import unittest

from denovonear.weights import get_distances, geomean, WeightedChoice, \
    analyse_de_novos, simulate_distribution

class TestSimulationsPy(unittest.TestCase):
    """ unit test the simulation functions
    """
    
    def setUp(self):
        """
        """
        
        # set up a range of possible positions, all with a uniform probability
        # of being selected
        self.choices = WeightedChoice()
        for x in range(1000):
            self.choices.add_choice(x, 0.0001)
        
        self.iterations = 100000
    
    def test_analyse_de_novos_dispersed(self):
        """ test analyse_de_novos() works correctly for dispersed de novos
        """
        
        # spread sites throughout a 1000 bp transcript
        positions = [100, 300, 600]
        distances = get_distances(positions)
        observed = geomean(distances)
        
        p_val = analyse_de_novos(self.choices, self.iterations, len(positions), observed)
        
        self.assertAlmostEqual(p_val, 0.635, places=2)
    
    def test_analyse_de_novos_clustered(self):
        """ test analyse_de_novos() works correctly for clustered de novos
        """
        
        # cluster sites within 20 bp in a 1000 bp transcript
        positions = [100, 110, 120]
        distances = get_distances(positions)
        observed = geomean(distances)
        
        p_val = analyse_de_novos(self.choices, 1000000, len(positions), observed)
        
        self.assertAlmostEqual(p_val, 0.002, places=3)
    
    def test_simulate_distribution(self):
        ''' check that simulate_distribution works correctly
        '''
        
        # repeated function calls should give different samples
        first = simulate_distribution(self.choices, iterations=5, de_novos_count=3)
        second = simulate_distribution(self.choices, iterations=5, de_novos_count=3)
        
        self.assertNotEqual(first, second)
        
