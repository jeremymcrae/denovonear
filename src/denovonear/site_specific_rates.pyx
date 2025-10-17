# cython: language_level=3, boundscheck=False
""" get weight gene mutation rates
"""

from libc.stdint cimport uint8_t, uint64_t
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from libcpp.string cimport string, stod
from libcpp cimport bool
from cython.operator cimport dereference as deref

from denovonear.weights cimport WeightedChoice
from gencodegenes.transcript cimport Region, Transcript, Codon

cdef unordered_map[string, unordered_map[string, double]] prepare_context_rates(rates: vector[vector[string]]):
    ''' convert list of list of strings [['AAA', 'AGA', '5e-8'], ...] to cpp object
    '''
    cdef  unordered_map[string, unordered_map[string, double]] parsed
    ref: string
    alt: string
    for row in rates:
        ref = row[0]
        alt = row[1]
        parsed[ref][alt] = stod(row[2])
    
    return parsed

cdef unordered_map[uint64_t, unordered_map[char, double]] prepare_site_rates(rates: dict[int, dict[str, double]]):
    ''' convert dict of {pos: {allele1: rate, ...}, ...} to matching cpp object
    '''
    cdef unordered_map[uint64_t, unordered_map[char, double]] parsed
    cdef unordered_map[char, double] site
    for pos, data in rates.items():
        site.clear()
        for allele, rate in data.items():
            site[ord(allele[0])] = rate
        parsed[pos] = site
    return parsed

cdef Region get_gene_range_cpp(tx: Transcript):
    ''' get the lowest and highest positions of a transcripts coding sequence
    '''
    start: int = min(tx.cds_start, tx.cds_end)
    end: int = max(tx.cds_start, tx.cds_end)
    
    return Region(start, end)

def get_gene_range(tx: Transcript):
    ''' get the lowest and highest positions of a transcripts coding sequence
    
    wrapper to cpp function
    '''
    region = get_gene_range_cpp(tx)
    return {'start': region.start, 'end': region.end}

cdef char get_mutated_aa_cpp(tx: Transcript,
                             codon: string,
                             char base,
                             intra_codon: int,
                             ):
    ''' find the amino acid resulting from a base change to a codon
    
    Args:
        tx: transcript object for a gene
        codon: DNA sequence of a single codon
        base: alternate base (e.g. 'G') to introduce
        intra_codon: position within the codon to be altered (0-based)
    
    Returns:
        single character amino acid code translated from the altered codon.
    '''
    # figure out what the mutated codon is
    codon[intra_codon] = base
    return tx.translate(codon.decode()).encode()[0]

def get_mutated_aa(tx: Transcript,
                   codon: str,
                   base: str,
                   intra_codon: int,
                   ):
    ''' find the amino acid resulting from a base change to a codon
    
    wrapper to cpp function
    '''
    return chr(get_mutated_aa_cpp(tx, codon.encode(), ord(base), intra_codon))

cdef str check_consequence(char initial_aa,
                           char mutated_aa,
                           offset: int,
                           boundary_dist: int
                           ):
    ''' get the consequence of an amino acid change (or not)
    '''
    
    cdef char stop = b'*'
    in_coding: bool = offset == 0
    if (initial_aa != stop) and (mutated_aa == stop):
        # checks if two amino acids are a nonsense (eg stop_gained) mutation
        return "nonsense"
    elif not in_coding and (boundary_dist < 3):
        # check if a variant has a splice_donor or splice_acceptor consequence
        # These variants are defined as being the two intronic base-pairs
        # adjacent to the intron/exon boundary.
        return "splice_lof"
    elif initial_aa != mutated_aa:
        # include the site if it mutates to a different amino acid.
        return "missense"
    elif in_coding:
        if boundary_dist < 4:
            # catch splice region variants within the exon, and in the appropriate
            # region of the intron (note that loss of function splice_donor and
            # splice_acceptor variants have been excluded when we spotted nonsense).
            return "splice_region"
    else:
        if boundary_dist < 9:
            # check for splice_region_variant inside intron
            return "splice_region"
        else:
            return "intronic"
    
    return "synonymous"

cdef class SiteRates:
    ''' class to run through all mutations for a transcript CDS, and figure out
    mutation rate and type. Allocates to a random sampler for the correct
    consequence type.
    '''
    cdef Region gene_range
    cdef unordered_map[char, char] transdict
    cdef bool has_mask
    cdef vector[char] bases
    cdef Transcript tx
    cdef Transcript masked
    
    cdef unordered_map[string, unordered_map[string, double]] context_rates
    cdef unordered_map[uint64_t, unordered_map[char, double]] per_pos_rates
    
    cdef int kmer_length
    cdef int mid_pos
    cdef object _choices
    cdef bool context_based
    cdef bool use_cds_coords
    
    def __cinit__(self,
                  transcript: Transcript,
                  rates,
                  masked_sites: Transcript | None=None,
                  cds_coords: bool=True,
                  ):
        '''
        '''
        
        self.context_based = isinstance(rates, list)
        self.use_cds_coords = cds_coords
        
        if self.context_based:
            self.context_rates = prepare_context_rates(rates)
        elif isinstance(rates, dict):
            self.per_pos_rates = prepare_site_rates(rates)
        else:
            raise ValueError('unknown format for mutation rates')
        
        if transcript is None:
            raise ValueError('no transcript supplied')
        
        cdef char A = b'A'
        cdef char C = b'C'
        cdef char G = b'G'
        cdef char T = b'T'
        self.transdict =  {A: T, T: A, G: C, C: G}
        self.bases = [A, C, G, T]
        
        self.tx = transcript
        self.has_mask = False
        if masked_sites is not None:
            self.masked = masked_sites
            self.has_mask = True
        
        self.clear()
        self.gene_range = get_gene_range_cpp(self.tx)
        self.process_cds()
    
    def clear(self):
        ''' construct empty WeightedChoice objects for each consequence type
        '''
        self._choices = {}
        for cq in ["missense", "nonsense", "synonymous", "splice_lof",
                   "splice_region", "loss_of_function", "intronic"]:
            self._choices[cq] = WeightedChoice()
    
    def process_cds(self):
        ''' run through all positions in the CDS and add SNVs to appropriate rates
        '''
        self.kmer_length = 3
        self.mid_pos = 1
        if self.context_based:
            # calculate the length of sequences used in the mutation rate dictionary
            # This means we can flexibly use 3-mers, 5-mers, 7-mers etc if desired.
            seq: string = deref(self.context_rates.begin()).first
            self.kmer_length = seq.length()
            self.mid_pos = self.kmer_length // 2
        
        # check the consequence alternates for each base in the coding sequence
        for i in range(self.gene_range.start, self.gene_range.end + 1):
            self.check_position(i)
    
    def __getitem__(self, category: str):
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
        return self._choices[category]
    
    cpdef check_position(self, bp: int):
        ''' add the consequence specific rates for the alternates for a variant
        
        Args:
            bp: genomic position of the variant
        '''
        # ignore sites within masked regions (typically masked because the
        # site has been picked up on alternative transcript)
        if self.has_mask and self.masked.in_coding_region(bp):
            return 
        
        #  ignore sites outside the CDS region
        if bp < self.gene_range.start or bp > self.gene_range.end:
            return
        
        boundary_dist: int = self.tx.get_boundary_distance(bp)
        seq: string = self.tx.get_centered_sequence(bp, self.kmer_length).encode()
        
        if self.tx.strand != '+':
            seq = self.tx.reverse_complement(seq.decode()).encode()
    
        try:
            codon = self.tx.get_codon_info(bp)
        except ValueError:
            return
        
        cdef char initial_aa = b'0'
        codon_seq: string= b''
        intra_codon: int= -1
        if codon['initial_aa'] is not None:
            initial_aa = codon['initial_aa'].encode()[0]
            codon_seq: string = codon['codon_seq'].encode()
            intra_codon: int = codon['intra_codon']
        
        codon_offset: int = codon['offset']
        cds_pos: int = codon['cds_pos']
    
        # drop the initial base, since we want to mutate to other bases
        cdef char ref = seq[self.mid_pos]
        cdef vector[char] alts = [x for x in self.bases if x != ref]
        
        if self.tx.strand != '+':
            ref = self.transdict[ref]
        
        cdef double rate
        cdef char mutated_aa
        alt_seq: string = seq
        category: str
        for alt in alts:
            mutated_aa = initial_aa
            alt_seq[self.mid_pos] = alt
            
            if initial_aa != b'0':
                mutated_aa = get_mutated_aa_cpp(self.tx, codon_seq, alt, intra_codon)
            
            category = check_consequence(initial_aa, mutated_aa, codon_offset, boundary_dist)
            
            # figure out the ref and alt alleles, with respect to the + strand.
            if self.tx.strand != '+':
                alt = self.transdict[alt]

            if self.context_based:
                rate = self.context_rates[seq][alt_seq]
            else:
                rate = self.per_pos_rates[bp][alt]
            
            if self.use_cds_coords:
                self._choices[category].add_choice(cds_pos, rate, chr(ref), chr(alt), codon_offset)
            else:
                self._choices[category].add_choice(bp, rate, chr(ref), chr(alt), 0)
            
            if (category == "nonsense") or (category == "splice_lof"):
                self._choices["loss_of_function"].add_choice(cds_pos, rate, chr(ref), chr(alt), codon_offset)
    
    cpdef check_consequence(self, initial_aa: str, mutated_aa: str, position: int):
        ''' check the consequence of an aminoa cid change at a given position
        
        This is just for testing purposes, so we can validate the actual
        check_consequence function
        '''
        if len(initial_aa) == 0:
            initial_aa = '0'
        if len(mutated_aa) == 0:
            mutated_aa = '0'
        
        cdef char initial_aa_char = initial_aa.encode()[0]
        cdef char mutated_aa_char = mutated_aa.encode()[0]
        boundary_dist = self.tx.get_boundary_distance(position)
        codon = self.tx.get_codon_info(position)
        return check_consequence(initial_aa_char, mutated_aa_char, codon['offset'], boundary_dist)
