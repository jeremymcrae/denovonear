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

from libcpp.vector cimport vector
from libcpp cimport bool
from cython.operator cimport dereference as deref

cdef extern from "weighted_choice.h":
    cdef cppclass Chooser:
        Chooser() except +
        void add_choice(int, double)
        int choice()
        double get_summed_rate()

cdef class WeightedChoice:
    cpdef Chooser *thisptr  # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.thisptr = new Chooser()
    def __dealloc__(self):
        del self.thisptr
    def add_choice(self, site, prob):
        """ add another possible choice for selection
        
        Args:
            site: an CDS position so we know each choice.
            prob: probability of selecting this choice.
        """
        self.thisptr.add_choice(site, prob)
    
    def choice(self):
        """ chooses a random element using a set of probability weights
        
        Returns:
            the name of the randomly selected element (e.g. position)
        """
        
        return self.thisptr.choice()
    
    def get_summed_rate(self):
        """ return the cumulative probability for the object
        """
        return self.thisptr.get_summed_rate()

cdef extern from "simulate.h":
    vector[double] _get_distances(vector[int])
    bool _has_zero(vector[double])
    double _geomean(vector[double])
    bool _halt_permutation(double, int, double, double)
    # vector[double] _simulate_distribution(Chooser, int, int)
    double _analyse_de_novos(Chooser, int, int, double)

def get_distances(vector[int] positions):
    """ gets the distances between two or more CDS positions
    
    Args:
        positions: list of CDS positions as ints
    
    Returns:
        list of pairwise distances between sites
    """
    
    return _get_distances(positions)

def has_zero(vector[double] distances):
    """ figure out whether any of the pairwise distances is zero
    """
    
    return _has_zero(distances)

def geomean(vector[double] distances):
    """ gets the geometric mean distance between two or more CDS positions
    
    Args:
        distances: list of distances between CDS positions
    
    Returns:
        provides the mean distance of the pairwise distances
    """
    
    return _geomean(distances)

# def simulate_distribution(WeightedChoice choices, int iterations, int de_novos_count):
#     """
#     """
#
#     return _simulate_distribution(deref(choices.thisptr), iterations, de_novos_count)
#
def analyse_de_novos(WeightedChoice choices, int iterations, int de_novos_count, double observed_value):
    """
    """

    return _analyse_de_novos(deref(choices.thisptr), iterations, de_novos_count, observed_value)
