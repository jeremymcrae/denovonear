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

from denovonear.weights import (get_distances,
                                get_structure_distances,
                                geomean,
                                geomean_double,
                                WeightedChoice,
                                simulate_distances,
                                simulate_structure_distances,
                                simulate_clustering,
                                simulate_structure_clustering,
                                )

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
    
    def test_simulate_clustering_dispersed(self):
        """ test simulate_clustering() works correctly for dispersed de novos
        """
        
        # spread sites throughout a 1000 bp transcript
        positions = [100, 300, 600]
        distances = get_distances(positions)
        observed = geomean(distances)
        
        p_val = simulate_clustering(self.choices, self.iterations, len(positions), observed)
        
        self.assertAlmostEqual(p_val, 0.635, delta=0.01)
    
    def test_simulate_clustering_clustered(self):
        """ test simulate_clustering() works correctly for clustered de novos
        """
        
        # cluster sites within 20 bp in a 1000 bp transcript
        positions = [100, 110, 120]
        distances = get_distances(positions)
        observed = geomean(distances)
        
        p_val = simulate_clustering(self.choices, 1000000, len(positions), observed)
        
        self.assertAlmostEqual(p_val, 0.002, places=3)
    
    def test_simulate_distances(self):
        ''' check that simulate_distances works correctly
        '''
        
        # repeated function calls should give different samples
        first = simulate_distances(self.choices, iterations=5, de_novos_count=3)
        second = simulate_distances(self.choices, iterations=5, de_novos_count=3)
        
        self.assertNotEqual(first, second)
        
    def test_simulate_structure_distances(self):
        ''' check that simulate_structure_distances works correctly
        '''
        
        coords = [{'x': float(i), 'y': float(i), 'z': float(i)} for i in range(len(self.choices) // 3)]
        
        # repeated function calls should give different samples
        first = simulate_structure_distances(self.choices, coords, iterations=5, de_novos_count=3)
        second = simulate_structure_distances(self.choices, coords, iterations=5, de_novos_count=3)
        
        self.assertNotEqual(first, second)
    
    def test_simulate_structure_clustering_dispersed(self):
        """ test simulate_structure_clustering() works correctly for dispersed de novos
        """
        
        # spread sites throughout a 1000 bp transcript
        dnms = [{'x': 30.0, 'y': 30.0, 'z': 30.0},
                {'x': 100.0, 'y': 100.0, 'z': 100.0},
                {'x': 200.0, 'y': 200.0, 'z': 200.0},
                ]
        distances = get_structure_distances(dnms)
        observed = geomean_double(distances)
        
        coords = [{'x': float(i), 'y': float(i), 'z': float(i)} for i in range(len(self.choices) // 3)]
        p_val = simulate_structure_clustering(
            self.choices, coords, self.iterations, len(dnms), observed)
        
        self.assertAlmostEqual(p_val, 0.649, delta=0.01)
    
    def test_simulate_structure_clustering_clustered(self):
        """ test simulate_structure_clustering() works correctly for clustered de novos
        """
        
        # cluster sites within 20 bp in a 1000 bp transcript
        dnms = [{'x': 101.0, 'y': 101.0, 'z': 101.0},
                {'x': 110.0, 'y': 110.0, 'z': 110.0},
                {'x': 111.0, 'y': 111.0, 'z': 111.0},
                ]
        distances = get_structure_distances(dnms)
        observed = geomean_double(distances)
        
        coords = [{'x': i, 'y': i, 'z': i} for i in range(len(self.choices) // 3)]
        p_val = simulate_structure_clustering(
            self.choices, coords, self.iterations, len(dnms), observed)
        
        self.assertAlmostEqual(p_val, 0.002, delta=0.001)
    
