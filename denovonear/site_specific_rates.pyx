""" get weight gene mutation rates
"""

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
from cython.operator cimport dereference as deref

from denovonear.weights cimport Chooser, WeightedChoice
from denovonear.transcript cimport Tx, Transcript, Region

cdef extern from "site_rates.h":
    cdef cppclass SitesChecks:
        SitesChecks(Tx, vector[vector[string]]) except +
        
        void initialise_choices()
        Chooser * __getitem__(string) except +
        
        bool splice_lof_check(string, string, int)
        bool nonsense_check(string, string, int)
        bool missense_check(string, string, int)
        bool splice_region_check(string, string, int)
        bool synonymous_check(string, string, int)
        
        void check_position(int bp)
    
    cdef Region _get_gene_range(Tx)
    cdef string _get_mutated_aa(Tx, string, string, int) except +

cdef class SiteRates:
    cdef SitesChecks *_checks  # hold a C++ instance which we're wrapping
    cdef dict categories
    def __cinit__(self, Transcript transcript, vector[vector[string]] rates):
        self._checks = new SitesChecks(deref(transcript.thisptr), rates)
        self.categories = {}
        
        # self._checks.set_mask(deref(masked_sites.thisptr))
        # self._checks = new SitesChecks(deref(transcript.thisptr), rates)
        # if masked is None:
        #     self._checks = new SitesChecks(deref(transcript.thisptr), rates)
        # else:
        #     self._checks = new SitesChecks(deref(transcript.thisptr), rates, deref(masked_sites.thisptr))
    
    def __dealloc__(self):
        del self._checks
        for x in self.categories:
            self.categories[x].__dealloc__()
            del self.categories[x]
    
    def __getitem__(self, category):
        ''' get site-specific mutation rates for each CDS base
    
        Args:
            category: string to indicate the consequence type. The permitted
                types are "missense", "nonsense", "synonymous",
                "loss_of_function", "splice_lof", and "splice_region".
        
        Returns:
            A WeightedChoice object for the CDS, where each position is paired
            with its mutation rate. We can then randomly sample sites from the
            CDS WeightedChoice object according to the probability of each site
            being mutated to the specific consequence type.
        '''
        
        cdef Chooser * chooser
        
        if category not in self.categories:
            chooser = self._checks.__getitem__(category)
            choices = WeightedChoice()
            choices.thisptr = chooser
            self.categories[category] = choices
        
        return self.categories[category]
    
    def clear(self):
        self._checks.initialise_choices()
    
    def splice_lof_check(self, initial_aa, mutated_aa, position):
        return self._checks.splice_lof_check(initial_aa, mutated_aa, position)
    
    def nonsense_check(self, initial_aa, mutated_aa, position):
        return self._checks.nonsense_check(initial_aa, mutated_aa, position)
    
    def missense_check(self, initial_aa, mutated_aa, position):
        return self._checks.missense_check(initial_aa, mutated_aa, position)
        
    def splice_region_check(self, initial_aa, mutated_aa, position):
        return self._checks.splice_region_check(initial_aa, mutated_aa, position)
        
    def synonymous_check(self, initial_aa, mutated_aa, position):
        return self._checks.synonymous_check(initial_aa, mutated_aa, position)
    
    def check_position(self, bp):
        self._checks.check_position(bp)

def get_gene_range(Transcript tx):
    region = _get_gene_range(deref(tx.thisptr))
    return {"start": region.start, "end": region.end}

def get_mutated_aa(Transcript tx,  base, string codon, int intra_codon):
    return _get_mutated_aa(deref(tx.thisptr), base, codon, intra_codon)
