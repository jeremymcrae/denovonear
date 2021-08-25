# cython: language_level=3, boundscheck=False, emit_linenums=True

import bisect
import logging

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
from libcpp.map cimport map

from pyfaidx import Fasta

from denovonear.transcript cimport Tx, Region, CDS_coords
from denovonear.transcript import Transcript

cdef extern from "gtf.h" namespace "gencode":
    cdef struct GTFLine:
        string chrom
        string feature
        int start
        int end
        string strand
        string symbol
        string tx_id
        string transcript_type
        bool is_principal
        
    GTFLine parse_gtfline(string line)

cdef extern from "gencode.h" namespace "gencode":
    cdef struct NamedTx:
        string symbol
        Tx tx
        bool is_principal
    
    vector[NamedTx] open_gencode(string, bool)

cpdef _parse_gtfline(string line):
    ''' python function for unit testing GTF parsing
    '''
    return parse_gtfline(line)

cdef _convert_exons(vector[Region] exons):
    ''' convert vector of exon Regions to list of lists
    
    We need exons and CDS as lists of lists for constructing the python 
    Transcript object.
    '''
    return [[y.start, y.end] for y in exons]

cpdef _open_gencode(gtf_path, coding_only=True):
    ''' python function for unit testing loading transcripts from GTF
    '''
    cdef vector[NamedTx] _transcripts = open_gencode(gtf_path.encode('utf8'), coding_only)
    
    transcripts = []
    for x in _transcripts:
        tx = x.tx
        chrom = tx.get_chrom().decode('utf8')
        start = tx.get_start()
        end = tx.get_end()
        exons = _convert_exons(tx.get_exons())
        cds = _convert_exons(tx.get_cds())
        strand = chr(tx.get_strand())
        tx_id = tx.get_name().decode('utf8')
        transcript = Transcript(tx_id, chrom, start, end, strand, exons, cds, offset=0)
        transcripts.append((x.symbol.decode('utf8'), transcript, x.is_principal))
    return transcripts

__genome_ = None

cdef class Gene:
    cdef string _symbol
    cdef vector[Tx] _transcripts
    cdef vector[bool] _principal
    cdef str _chrom
    cdef int _start, _end
    def __cinit__(self, string symbol):
        self._symbol = symbol
        self.start = 999999999
        self.end = -999999999
    
    cdef add_tx(self, Tx tx, bool is_principal):
        self._transcripts.push_back(tx)
        self._principal.push_back(is_principal)
        self.chrom = tx.get_chrom().decode('utf8')
        self.start = min(self.start, tx.get_start())
        self.end = max(self.end, tx.get_end())
    
    def __repr__(self):
        chrom = self.chrom
        return f'Gene("{self.symbol}", {chrom}:{self.start}-{self.end})'
    
    @property
    def symbol(self):
        return self._symbol.decode('utf8')
    
    @property
    def chrom(self):
        return self._chrom
    @chrom.setter
    def chrom(self, value):
        self._chrom = value
    
    @property
    def start(self):
        return self._start
    @start.setter
    def start(self, value):
        self._start = value
    
    @property
    def end(self):
        return self._end
    @end.setter
    def end(self, value):
        self._end = value
    
    @property
    def strand(self):
        if self._transcripts.size() > 0:
            return chr(self._transcripts[0].get_strand())
        raise IndexError('no transcripts in gene yet')
    
    cdef _convert_exons(self, vector[Region] exons):
        ''' convert vector of exon Regions to list of lists
        
        We need exons and CDS as lists of lists for constructing the python 
        Transcript object.
        '''
        return [[y.start, y.end] for y in exons]
    
    cdef _to_Transcript(self, Tx tx):
        ''' construct Transcript (python object) from Tx (c++ object)
        '''
        offset = 10
        start = tx.get_start()
        end = tx.get_end()
        exons = self._convert_exons(tx.get_exons())
        cds = self._convert_exons(tx.get_cds())
        seq = __genome_[self.chrom][start-1-offset:end-1+offset].seq
        strand = self.strand
        tx_id = tx.get_name().decode('utf8')
        return Transcript(tx_id, self.chrom, start, end, strand, exons, cds, seq, offset=offset)
    
    @property
    def transcripts(self):
        ''' get list of Transcripts for gene, with genomic DNA included
        '''
        return [self._to_Transcript(x) for x in self._transcripts]
    
    cdef int _cds_len(self, Tx tx):
        ''' get length of coding sequence for a Tx object based transcript
        '''
        cdef CDS_coords coords = tx.get_coding_distance(tx.get_cds_end())
        return coords.position + 1
    
    cdef Tx _max_by_cds(self, vector[Tx] transcripts):
        ''' get longest transcript by CDS length
        '''
        cdef Tx max_tx
        length = 0
        for tx in transcripts:
            curr_len = self._cds_len(tx)
            if curr_len > length:
                length = curr_len
                max_tx = tx
        return max_tx
    
    @property
    def canonical(self):
        ''' find the canonical transcript for a gene.
        
        Canonical is defined as:
            - transcript with longest CDS tagged with appris_principal in the GTF
            - if no appris_principal tags for any tx, use the tx with longest CDS
        
        Occasionally there are multiple transcripts tagged as appris_principal
        and with the same longest CDS, we use the first one of those.
        '''
        cdef vector[Tx] principal
        for i in range(self._transcripts.size()):
            if self._principal[i]:
                principal.push_back(self._transcripts[i])
        
        if principal.size() == 0:
            principal = self._transcripts
        
        # cdef Tx tx = self._max_by_cds(principal)
        return self._to_Transcript(self._max_by_cds(principal))
    
    def in_any_tx_cds(self, pos):
        ''' find if a pos is in coding region of any transcript of a gene
        '''
        return any(tx.in_coding_region(pos) for tx in self._transcripts)

cdef class Gencode:
    cdef dict genes
    cdef list coords
    def __cinit__(self, gencode, fasta, coding_only=True):
        ''' initialise Gencode
        
        Args:
            gencode: path to gencode annotations file
            fasta: path to fasta for genome matching annotations build
            coding: restrict to protein_coding only by default
        '''
        logging.info(f'opening genome fasta: {fasta}')
        global __genome_
        __genome_ = Fasta(str(fasta))
        logging.info(f'opening gencode annotations: {gencode}')
        cdef vector[NamedTx] transcripts = open_gencode(str(gencode).encode('utf8'), coding_only)
        cdef Gene curr
        self.genes = {}
        for x in transcripts:
            symbol = x.symbol.decode('utf8')
            if symbol not in self.genes:
                self.genes[symbol] = Gene(symbol.encode('utf8'))
            curr = self.genes[symbol]
            curr.add_tx(x.tx, x.is_principal)
            self.genes[symbol] = curr
        
        self.coords = sorted(((x.chrom, x.start, x.end), symbol) for symbol, x in self.genes.items())
    
    def __repr__(self):
        return f'Gencode(n_genes={len(self)})'
    def __len__(self):
        return len(self.genes)
    def __getitem__(self, symbol):
        return self.genes[symbol]
    
    def distance(self, gene, chrom, pos):
        ''' get distance to nearest boundary of a gene
        '''
        if gene.chrom != chrom:
            return None
        if gene.start <= pos <= gene.end:
            return 0
        return min(abs(gene.start - pos), abs(gene.end - pos))
    
    def nearest(self, chrom, pos):
        ''' find the nearest gene to a genomic chrom, pos coordinate
        '''
        chrom = f'chr{chrom}' if not chrom.startswith('chr') else chrom
        i = bisect.bisect_left(self.coords, ((chrom, pos, pos), 'AAAA'))
        i = min(i, len(self.coords) - 1)  # ensure we don't go past the last gene
        
        if chrom != self[self.coords[i][1]].chrom:
            # the gene coords are in a single sorted list (by chrom and pos), and
            # we can match to the next chrom, just step back to the previous chrom
            i -= 1
        
        if pos < self[self.coords[i][1]].start:
            left = self[self.coords[max(i - 1, 0)][1]]
            right = self[self.coords[i][1]]
        else:
            left = self[self.coords[i][1]]
            right = self[self.coords[min(i + 1, len(self.coords) - 1)][1]]
        
        left_delta = self.distance(left, chrom, pos)
        right_delta = self.distance(right, chrom, pos)
        
        if right_delta is None and left_delta is None:
            raise ValueError(f"can't find any genes on {chrom}")
        if right_delta is None:
            return left
        elif left_delta is None:
            return right
        elif left_delta < right_delta:
            return left
        else:
            return right
    
    def in_region(self, chrom, start, end):
        ''' find genes within a genomic region
        
        Args:
            chrom: chromosome to search on
            start: start position of region
            end: end position of region
        
        Returns:
            list of Gene objects
        '''
        chrom = f'chr{chrom}' if not chrom.startswith('chr') else chrom
        left = bisect.bisect_left(self.coords, ((chrom, start, start), 'AAAA'))
        right = bisect.bisect_left(self.coords, ((chrom, end, end), 'AAAA'))
        
        genes = (self.genes[self.coords[i][1]] for i in range(left, right))
        return [x for x in genes if x.chrom == chrom]
    
    def __exit__(self):
        ''' cleanup at exit
        '''
        global __genome_
        __genome_.close()
        __genome_ = None