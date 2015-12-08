# distutils: language = c++
# distutils: sources = weighted_choice.cpp

cdef extern from "weighted_choice.h" namespace "weights":
    cdef cppclass WeightedChoice:
        WeightedChoice() except +
        int add_choice(int, double)
        int choice()
        int get_summed_rate()

cdef class PyWeightedChoice:
    cdef WeightedChoice *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.thisptr = new WeightedChoice()
    def __dealloc__(self):
        del self.thisptr
    def add_choice(self, site, prob):
        return self.thisptr.add_choice(site, prob)
    def choice(self):
        return self.thisptr.choice()
    def get_summed_rate(self):
        self.thisptr.get_summed_rate()
