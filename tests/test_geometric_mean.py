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

from denovonear.weights import get_distances, geomean

class TestGeomeanPy(unittest.TestCase):
    """ unit test the geomean function
    """
    
    def test_geomean(self):
        """ test geomean() works correctly
        """
        
        self.assertEqual(geomean([0]), 0)
        self.assertEqual(geomean([1]), 1)
        self.assertEqual(geomean([1, 1]), 1)
        self.assertEqual(geomean([1, 2]), 1.4142135623730951)
        self.assertEqual(geomean([10]), 10)
        self.assertEqual(geomean([0, 900, 900]), 92.2860120092046)
        self.assertEqual(geomean([10, 890, 900]), 200.08329863520368)
        
        # check that if we try to get the mean distance for a list with fewer
        # than two positions in it, we get an error, since we can't estimate any
        # distances
        self.assertTrue(math.isnan(geomean([])))
    
    def test_get_distance(self):
        """ test that get_distances() works correctly
        """
        
        # check pairwise distances
        self.assertEqual(get_distances([0, 0]), [0])
        self.assertEqual(get_distances([0, 1]), [1])
        self.assertEqual(get_distances([0, 1, 10]), [1, 10, 9])
        
        # empty lists (or lists with one entry) can't give pairwise distances
        self.assertEqual(get_distances([]), [])
        self.assertEqual(get_distances([0]), [])
        
        # check that float positionss are automatically converted to integers
        self.assertEqual(get_distances([0.5, 10]), [10])
        
        # check that negative positions still work
        self.assertEqual(get_distances([-10, -5]), [5])
        
        # and check we get an error if we call with strings
        with self.assertRaises(TypeError):
            get_distances([0, 1, "e"])
        
    
