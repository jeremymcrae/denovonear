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
import itertools

def geomean(positions):
    """ gets the geometric mean distance between two or more CDS positions
    
    Args:
        positions: list of numbers
    
    Returns:
        provides the mean distance of the position pair distances
    """
    
    if len(positions) == 2:
        return abs(positions[0] - positions[1])
    
    # figure the distances between every possible pair
    distances = [ abs(x - y) for x, y in itertools.combinations(positions, 2) ]
    
    # get the geometric mean, but be careful around values of 0, since
    # without correction, the mean distance would be zero
    if 0 in distances:
        # allow for 0s in a geometric mean by shifting everything up one,
        # then dropping the mean by one at the end
        distances = [ math.log10(x + 1) for x in distances ]
        logmean = sum(distances)/len(distances)
        mean = (10 ** logmean) - 1
    else:
        distances = [ math.log10(x) for x in distances ]
        logmean = sum(distances)/len(distances)
        mean = 10 ** logmean
    
    return mean
