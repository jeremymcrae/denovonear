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

cdef class WeightedChoice:
    def __cinit__(self):
        self.thisptr = new Chooser()
        self.pos = 0
    
    def __dealloc__(self):
        if self.thisptr is not NULL:
            del self.thisptr
    
    def __len__(self):
        return self.thisptr.len()
    
    def __iter__(self):
        return self
    
    def __next__(self):
        if self.pos >= len(self):
            self.pos = 0
            raise StopIteration
        
        iter = self.thisptr.iter(self.pos)
        self.pos += 1
        
        return {"pos": iter.pos, "ref": iter.ref.decode('utf8'),
            "alt": iter.alt.decode('utf8'), 'prob': iter.prob, "offset": iter.offset}
    
    def append(self, WeightedChoice other):
        ''' combines the sites from two WeightedChoice objects
    
        NOTE: This is inefficient, but I can't figure out an easier way around
        NOTE: cython.
    
        Args:
            other: WeightedChoice object
        '''
    
        self.thisptr.append(deref(other.thisptr))
    
    def add_choice(self, site, prob, ref='N', alt='N', offset=0):
        """ add another possible choice for selection
        
        Args:
            site: a CDS position of the selected base.
            prob: probability of selecting this base.
            ref: string for reference allele
            alt: string for alternate allele
        """
        
        ref = ref.encode('utf8')
        alt = alt.encode('utf8')
        
        if len(ref) > 1 or len(alt) > 1:
            raise TypeError("requires single base alleles: {}, {}".format(ref, alt))
        
        self.thisptr.add_choice(site, prob, ref, alt, offset)
    
    def choice(self):
        """ chooses a random element using a set of probability weights
        
        Returns:
            the name of the randomly selected element (e.g. position)
        """
        
        return self.choice_with_alleles()["pos"]
    
    def choice_with_alleles(self):
        """ chooses a random element, but include alleles in output
        
        Returns:
            a dictionary of the randomly selected element, with "pos", "ref"
            and "alt" entries.
        """
        
        choice = self.thisptr.choice()
        
        return {"pos": choice.pos, "ref": choice.ref.decode('utf8'),
            "alt": choice.alt.decode('utf8'), "offset": choice.offset}
    
    def get_summed_rate(self):
        """ return the cumulative probability for the object
        """
        return self.thisptr.get_summed_rate()

cdef extern from "simulate.h":
    vector[int] _get_distances(vector[int])
    bool _has_zero(vector[int])
    double _geomean(vector[int])
    bool _halt_permutation(double, int, double, double)
    vector[double] _simulate_distribution(Chooser, int, int)
    double _analyse_de_novos(Chooser, int, int, double)

def get_distances(vector[int] positions):
    """ gets the distances between two or more CDS positions
    
    Args:
        positions: list of CDS positions as ints
    
    Returns:
        list of pairwise distances between sites
    """
    
    return _get_distances(positions)

def has_zero(vector[int] distances):
    """ figure out whether any of the pairwise distances is zero
    """
    
    return _has_zero(distances)

def geomean(vector[int] distances):
    """ gets the geometric mean distance between two or more CDS positions
    
    Args:
        distances: list of distances between CDS positions
    
    Returns:
        provides the mean distance of the pairwise distances
    """
    
    return _geomean(distances)

def simulate_distribution(WeightedChoice choices, int iterations, int de_novos_count):
    """
    """
    
    return _simulate_distribution(deref(choices.thisptr), iterations, de_novos_count)

def analyse_de_novos(WeightedChoice choices, int iterations, int de_novos_count, double observed_value):
    """
    """
    
    return _analyse_de_novos(deref(choices.thisptr), iterations, de_novos_count, observed_value)
