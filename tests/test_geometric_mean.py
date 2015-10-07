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

import unittest

from denovonear.geometric_mean import geomean

class TestGeomeanPy(unittest.TestCase):
    """ unit test the geomean function
    """
    
    def test_geomean(self):
        """ test get_position_in_cds() works correctly
        """
        
        self.assertEqual(geomean([100, 100]), 0)
        self.assertEqual(geomean([100, 110]), 10)
        self.assertEqual(geomean([100, 100, 1000]), 92.2860120092046)
        self.assertEqual(geomean([100, 110, 1000]), 200.08329863520368)
        
        # check that if we try to get the mean distance for a list with fewer
        # than two positions in it, we get an error, since we can't estimate any
        # distances
        self.assertRaises(ZeroDivisionError, geomean, [100])
        self.assertRaises(ZeroDivisionError, geomean, [])
        
    
