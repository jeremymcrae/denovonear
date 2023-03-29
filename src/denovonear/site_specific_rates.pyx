# cython: language_level=3, boundscheck=False
""" get weight gene mutation rates
"""

from libc.stdint cimport uint64_t
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from libcpp.string cimport string
from libcpp cimport bool
from cython.operator cimport dereference as deref

from denovonear.weights cimport Chooser, WeightedChoice
from gencodegenes.transcript cimport Tx, Transcript, Region, Codon

cdef extern from "site_rates.h":
    cdef cppclass SitesChecks:
        SitesChecks(Tx, vector[vector[string]], bool) except +
        SitesChecks(Tx, unordered_map[uint64_t, unordered_map[char, double]], bool) except +
        
        Tx _tx
        void initialise_choices()
        Chooser * __getitem__(string) except +
        
        void add_mask(Tx mask)
        void check_position(int)
        void check_consequence(string, char, char, int)
    
    cdef Region _get_gene_range(Tx)
    cdef char _get_mutated_aa(Tx, char, string, int) except +

cdef vector[vector[string]] prepare_rates_vector(rates):
    ''' convert list of list of strings [['AAA', 'AGA', '5e-8'], ...] to cpp object
    '''
    cdef vector[vector[string]] rates_vector
    cdef vector[string] row
    for x in rates:
        row.clear()
        for y in x:
            row.push_back(y)
        rates_vector.push_back(row)
    return rates_vector

cdef unordered_map[uint64_t, unordered_map[char, double]] prepare_rates_dict(rates):
    ''' convert dict of {pos: {allele1: rate, ...}, ...} to matching cpp object
    '''
    cdef unordered_map[uint64_t, unordered_map[char, double]] rates_dict 
    cdef unordered_map[char, double] site
    for pos, data in rates.items():
        site.clear()
        for allele, rate in data.items():
            site[ord(allele[0])] = rate
        rates_dict[pos] = site
    return rates_dict

cdef class SiteRates:
    cdef SitesChecks *_checks  # hold a C++ instance which we're wrapping
    def __cinit__(self, Transcript transcript, rates,
            Transcript masked_sites=None, bool cds_coords=True):
        
        cdef bool uses_context = isinstance(rates, list)
        cdef vector[vector[string]] rates_vector 
        cdef unordered_map[uint64_t, unordered_map[char, double]] rates_dict 
        if isinstance(rates, list):
            rates_vector = prepare_rates_vector(rates)
        elif isinstance(rates, dict):
            rates_dict = prepare_rates_dict(rates)
        else:
            raise ValueError('unknown format for mutation rates')

        if transcript is None:
            raise ValueError('no transcript supplied')
        
        if uses_context:
            self._checks = new SitesChecks(deref(transcript.thisptr), rates_vector,
                    cds_coords)
        else:
            self._checks = new SitesChecks(deref(transcript.thisptr), rates_dict,
                    cds_coords)
        
        if masked_sites:
            self._checks.add_mask(deref(masked_sites.thisptr))
    
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
    
    cpdef check_consequence(self, initial_aa, mutated_aa, position):
        if len(initial_aa) == 0:
            initial_aa = '0'
        if len(mutated_aa) == 0:
            mutated_aa = '0'
        cdef string cq = b"synonymous"
        initial_aa = ord(initial_aa)
        mutated_aa = ord(mutated_aa)
        self.check_position(position)
        codon = self._checks._tx.get_codon_info(position)
        self._checks.check_consequence(cq, initial_aa, mutated_aa, codon.offset)
        return cq.decode('utf8')

def get_gene_range(Transcript tx):
    region = _get_gene_range(deref(tx.thisptr))
    return {"start": region.start, "end": region.end}

def get_mutated_aa(Transcript tx,  base, codon, intra_codon):
    
    codon = codon.encode('utf8')
    
    return chr(_get_mutated_aa(deref(tx.thisptr), ord(base), codon, intra_codon))
