# distutils: language = c++
# distutils: sources = denovonear/cpp/simulate.cpp

# from denovonear.weights import WeightedChoice
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "weighted_choice.h":
    cdef cppclass Chooser:
        Chooser() except +
        void add_choice(int, double)
        int choice()
        double get_summed_rate()

cdef extern from "simulate.h":
    vector[double] get_distances(vector[int])
    bool has_zero(vector[double])
    double get_geomean(vector[double])
    vector[double] simulate_distribution(Chooser, int, int)
    bool halt_permutation(double, int, double, double)
    double analyse_de_novos(Chooser, int, int, double)

def geomean(distances):
    get_geomean(distances)
