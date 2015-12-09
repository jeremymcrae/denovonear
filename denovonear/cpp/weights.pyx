# distutils: language = c++
# distutils: sources = denovonear/cpp/weighted_choice.cpp

cdef extern from "weighted_choice.h":
    cdef cppclass Chooser:
        Chooser() except +
        void add_choice(int, double)
        int choice()
        double get_summed_rate()

cdef class WeightedChoice:
    cdef Chooser *thisptr  # hold a C++ instance which we're wrapping
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
