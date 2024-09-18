# cython: language_level=3, boundscheck=False
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

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "weighted_choice.h":
    cdef cppclass Chooser:
        Chooser() except +
        void add_choice(int, double, char, char, int)
        int choice_pos_only()
        AlleleChoice choice()
        double get_summed_rate()
        int len()
        AlleleChoice iter(int)
        void append(Chooser)
    
    cdef struct AlleleChoice:
        int pos
        char ref
        char alt
        double prob
        int offset

cdef extern from "simulate.h":
    cdef struct Coord:
        double x
        double y
        double z
    
    void _get_distances(vector[int] &, vector[int] &)
    void _get_structure_distances(vector[Coord] &, vector[int] &)
    double _geomean(vector[int] &)
    bool _halt_permutation(double, int, double, double)
    vector[double] _simulate_distribution(Chooser, int, int)
    vector[double] _simulate_structure_distribution(Chooser, vector[Coord], int, int)
    double _analyse_de_novos(Chooser, int, int, double)
    double _analyse_structure_de_novos(Chooser, vector[Coord], int, int, double)

cdef class WeightedChoice:
    cdef int pos
    cdef Chooser *thisptr # hold a C++ instance which we're wrapping
