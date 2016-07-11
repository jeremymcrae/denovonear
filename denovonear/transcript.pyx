'''
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
'''

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "tx.h":
    cdef cppclass Tx:
        Tx(string, string, int, int, char) except +
        
        void set_exons(vector[vector[int]], vector[vector[int]]) except +
        void set_cds(vector[vector[int]]) except +
        Region fix_cds_boundary(int) except +
        
        vector[Region] get_exons()
        vector[Region] get_cds()
        string get_name()
        string get_chrom()
        int get_start()
        int get_end()
        char get_strand()
        int get_cds_start()
        int get_cds_end()
        
        bool in_exons(int)
        Region find_closest_exon(int)
        Region find_closest_exon(int, vector[Region])
        bool in_coding_region(int)
        int get_exon_containing_position(int, vector[Region]) except +
        int get_coding_distance(int, int) except +
        int chrom_pos_to_cds(int) except +
    
    cdef struct Region:
        int start
        int end

cdef class Transcript:
    cpdef Tx *thisptr  # hold a C++ instance which we're wrapping
    def __cinit__(self, transcript_id, chrom, start, end, strand):
        self.thisptr = new Tx(transcript_id, chrom, start, end, ord(strand))
    def __dealloc__(self):
        del self.thisptr
    
    def __repr__(self):
        return "{0}({1} {2}:{3}-{4})".format(self.__class__.__name__,
            self.thisptr.get_name(), self.thisptr.get_chrom(),
            self.thisptr.get_start(), self.thisptr.get_end())
    
    def __str__(self):
        return self.__repr__()
    
    def __hash__(self):
        return hash((self.thisptr.get_chrom(), self.thisptr.get_start(),
            self.thisptr.get_end()))
    
    def set_exons(self, exon_ranges, cds_ranges):
        ''' add exon ranges
        
        Args:
            exon_ranges: a CDS position of the selected base.
        '''
        
        self.thisptr.set_exons(exon_ranges, cds_ranges)
    
    
    def set_cds(self, cds_ranges):
        ''' set CDS ranges
        
        Args:
            cds_ranges: a CDS position of the selected base.
        '''
        
        self.thisptr.set_cds(cds_ranges)
        
    def get_exons(self):
        return self.thisptr.get_exons()
    def get_cds(self):
        return self.thisptr.get_cds()
    def get_name(self):
        return self.thisptr.get_name()
    def get_chrom(self):
        return self.thisptr.get_chrom()
    def get_start(self):
        return self.thisptr.get_start()
    def get_end(self):
        return self.thisptr.get_end()
    def get_strand(self):
        return chr(self.thisptr.get_strand())
    def get_cds_start(self):
        return self.thisptr.get_cds_start()
    def get_cds_end(self):
        return self.thisptr.get_cds_end()
    def fix_cds_boundary(self, pos):
        return self.thisptr.fix_cds_boundary(pos)
    
    def in_exons(self, position):
        ''' check if a site lies within the exon ranges
        
        Args:
            position: an integer-based chromosome position.
        '''
        
        return self.thisptr.in_exons(position)
    
    def find_closest_exon(self, position, exons=None):
        ''' finds the positions of the exon closest to a position
        
        Can optionally specify a different set of exon ranges to check from,
        such as the CDS positions.
        '''
        
        cdef vector[Region] temp
        cdef Region pair
        
        if exons is None:
            return self.thisptr.find_closest_exon(position)
        else:
            for (start, end) in exons:
                pair.start = start
                pair.end = end
                temp.push_back(pair)
            return self.thisptr.find_closest_exon(position, temp)
    
    def in_coding_region(self, position):
        return self.thisptr.in_coding_region(position)
    
    def get_exon_containing_position(self, position, exons):
        '''
        '''
        
        cdef vector[Region] regions
        cdef Region pair
        
        for (start, end) in exons:
            pair.start = start
            pair.end = end
            regions.push_back(pair)
        
        return self.thisptr.get_exon_containing_position(position, regions)
    
    def get_coding_distance(self, pos_1, pos_2):
        return self.thisptr.get_coding_distance(pos_1, pos_2)
        
    def chrom_pos_to_cds(self, pos_1):
        return self.thisptr.chrom_pos_to_cds(pos_1)
