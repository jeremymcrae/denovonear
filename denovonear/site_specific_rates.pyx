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
        SitesChecks(Tx, vector[vector[string]], bool) except +
        SitesChecks(Tx, vector[vector[string]], bool, Tx) except +
        
        void initialise_choices()
        Chooser * __getitem__(string) except +
        
        void check_position(int)
        string check_consequence(string, string, int)
    
    cdef Region _get_gene_range(Tx)
    cdef string _get_mutated_aa(Tx, string, string, int) except +

cdef class SiteRates:
    cdef SitesChecks *_checks  # hold a C++ instance which we're wrapping
    def __cinit__(self, Transcript transcript, vector[vector[string]] rates,
            Transcript masked_sites=None, cds_coords=True):
        
        if transcript is None:
            raise ValueError('no transcript supplied')
        
        if masked_sites is None:
            self._checks = new SitesChecks(deref(transcript.thisptr), rates,
                cds_coords)
        else:
            self._checks = new SitesChecks(deref(transcript.thisptr), rates,
                cds_coords, deref(masked_sites.thisptr))
    
    def __dealloc__(self):
        del self._checks
    
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
        
        cdef Chooser * chooser = self._checks.__getitem__(category.encode('utf8'))
        
        choices = WeightedChoice()
        choices.thisptr.append(deref(chooser))
        
        return choices
    
    def clear(self):
        self._checks.initialise_choices()
    
    def check_position(self, bp):
        self._checks.check_position(bp)
    
    def check_consequence(self, initial_aa, mutated_aa, position):
        initial_aa = initial_aa.encode('utf8')
        mutated_aa = mutated_aa.encode('utf8')
        return self._checks.check_consequence(initial_aa, mutated_aa, position).decode('utf8')

def get_gene_range(Transcript tx):
    region = _get_gene_range(deref(tx.thisptr))
    return {"start": region.start, "end": region.end}

def get_mutated_aa(Transcript tx,  base, codon, intra_codon):
    
    base = base.encode('utf8')
    codon = codon.encode('utf8')
    
    return _get_mutated_aa(deref(tx.thisptr), base, codon, intra_codon).decode('utf8')
