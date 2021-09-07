# cython: language_level=3, boundscheck=False
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
        Tx() except +
        
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
        
        bool is_exonic(int)
        int closest_exon_num(int)
        Region get_closest_exon(int)
        bool in_coding_region(int)
        CDS_coords to_closest_exon(int, Region)
        CDS_coords get_coding_distance(int) except +
        
        int get_position_on_chrom(int, int) except +
        int get_codon_number_for_cds_position(int)
        int get_position_within_codon(int)
        void add_cds_sequence(string)
        void add_genomic_sequence(string, int) except +
        string get_cds_sequence()
        string get_genomic_sequence()
        int get_genomic_offset()
        
        string reverse_complement(string)
        string get_centered_sequence(int, int) except +
        string get_codon_sequence(int) except +
        string get_seq_in_region(int, int) except +
        string translate(string) except +
        
        Codon get_codon_info(int) except +
        int get_boundary_distance(int) except +
        string consequence(int, int, string) except +
    
    cdef struct CDS_coords:
        int position
        int offset
    
    cdef struct Region:
        int start
        int end
    
    cdef struct Codon:
        int cds_pos
        string codon_seq
        int intra_codon
        int codon_number
        char initial_aa
        int offset

cdef class Transcript:
    cdef Tx *thisptr # hold a C++ instance which we're wrapping
