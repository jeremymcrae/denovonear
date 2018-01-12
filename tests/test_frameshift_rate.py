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
import tempfile

from denovonear.frameshift_rate import include_frameshift_rates

class TestIncludeFrameshiftRatesPy(unittest.TestCase):
    """ unit test the geomean function
    """
    
    def test_include_frameshift_rates(self):
        """ test include_frameshift_rates() works correctly
        """
        
        temp = tempfile.NamedTemporaryFile(mode='wt')
        temp.write('transcript_id\tchrom\tlength\tmissense_rate\tnonsense_rate'
            '\tsplice_lof_rate\tsplice_region_rate\tsynonymous_rate\n')
        temp.write('A\tchr1\t100\t-1\t-1\t-1\t-1\t-1\n')
        temp.write('B\tchr1\t1000\t-1\t-1\t-1\t-1\t-1\n')
        temp.write('C\tchr1\t1000\t-1\t-4\t-1\t-1\t-1\n')
        temp.flush()
        
        include_frameshift_rates(temp.name)
        temp.flush()
        
        with open(temp.name) as handle:
            values = [ line.strip().split('\t')[-1] for line in handle ]
        
        header, a, b, c = values
        
        self.assertEqual(header, 'frameshift_rate')
        self.assertAlmostEqual(float(a), -1.9240621930896515)
        self.assertAlmostEqual(float(b), -0.9240621930896515)
        self.assertAlmostEqual(float(c), -0.9240621930896515)
    
