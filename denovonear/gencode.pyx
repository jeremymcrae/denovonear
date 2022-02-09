# cython: language_level=3, boundscheck=False, emit_linenums=True

import bisect
import logging

from libcpp.algorithm cimport lower_bound, upper_bound
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
    
    cdef struct GenePoint:
        int pos
        string symbol
    
    vector[NamedTx] open_gencode(string, bool)
    bool CompFunc(const GenePoint &l, const GenePoint &r)
    vector[string] _in_region(string chrom, int start, int end, 
        map[string, vector[GenePoint]] & starts, map[string, vector[GenePoint]] & ends,
        int max_window) except+

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
    def __cinit__(self, symbol):
        if isinstance(symbol, str):
            symbol = symbol.encode('utf8')
        self._symbol = symbol
        self.start = 999999999
        self.end = -999999999
    
    cdef add_tx(self, Tx tx, bool is_principal):
        self._transcripts.push_back(tx)
        self._principal.push_back(is_principal)
        self.chrom = tx.get_chrom().decode('utf8')
        self.start = min(self.start, tx.get_start())
        self.end = max(self.end, tx.get_end())
    
    def add_transcript(self, _tx):
        ''' add a Transcript to the gene object
        
        This ends up coping the data from the Transcript object, rather than
        reusing the Tx object contained in the Trnascript, but it's not too much
        time wasted, so long as we don't do this millions of times.
        '''
        assert isinstance(_tx, Transcript)
        # construct a new Tx obect by copying out the relevant data
        cdef string tx_id = _tx.get_name().encode('utf8')
        cdef string chrom = _tx.get_chrom().encode('utf8')
        cdef int start = _tx.get_start()
        cdef int end = _tx.get_end()
        cdef char strand = ord(_tx.get_strand())
        cdef vector[vector[int]] exons
        cdef vector[vector[int]] cds
        cdef vector[int] exon
        for x in _tx.get_exons():
            exon = [x['start'], x['end']]
            exons.push_back(exon)
        for x in _tx.get_cds():
            exon = [x['start'], x['end']]
            cds.push_back(exon)
        
        cdef Tx tx = Tx(tx_id, chrom, start, end, strand)
        tx.set_exons(exons, cds)
        tx.set_cds(cds)
        
        cdef string seq = _tx.get_genomic_sequence().encode('utf8')
        cdef int offset = _tx.get_genomic_offset()
        if len(seq) > 0:
            if chr(strand) == '-':
                _tx.reverse_complement(seq)
                seq = _tx.reverse_complement(seq).encode('utf8')
            tx.add_genomic_sequence(seq, offset)
        self.add_tx(tx, False)
    
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
        offset = 5 if tx.get_genomic_offset() == 0 else tx.get_genomic_offset()
        chrom = tx.get_chrom().decode('utf8')
        start = tx.get_start()
        end = tx.get_end()
        exons = self._convert_exons(tx.get_exons())
        cds = self._convert_exons(tx.get_cds())
        seq = tx.get_genomic_sequence().decode('utf8')
        if seq == '':
            seq = None
        if seq is None and __genome_ is not None:
            seq = __genome_[chrom][start-1-offset:end-1+offset].seq.upper()
        strand = chr(tx.get_strand())
        if strand == '-' and seq is not None:
            seq = tx.reverse_complement(seq.encode('utf8')).decode('utf8')
        tx_id = tx.get_name().decode('utf8')
        return Transcript(tx_id, chrom, start, end, strand, exons, cds, seq, offset=offset)
    
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
    
    def distance(self, chrom, pos):
        ''' get distance to nearest boundary of a gene
        '''
        # sanatize the chromosome first
        if self.chrom.startswith('chr') and not chrom.startswith('chr'):
            chrom = f'chr{chrom}'
        elif not self.chrom.startswith('chr') and chrom.startswith('chr'):
            chrom = chrom[3:]
        
        if self.chrom != chrom:
            return None
        if self.start <= pos <= self.end:
            return 0
        return min(abs(self.start - pos), abs(self.end - pos))

cdef class Gencode:
    cdef dict genes
    cdef map[string, vector[GenePoint]] starts, ends
    def __cinit__(self, gencode=None, fasta=None, coding_only=True):
        ''' initialise Gencode
        
        Args:
            gencode: path to gencode annotations file
            fasta: path to fasta for genome matching annotations build
            coding: restrict to protein_coding only by default
        '''
        self.genes = {}
        if fasta:
            logging.info(f'opening genome fasta: {fasta}')
            global __genome_
            __genome_ = Fasta(str(fasta))
        logging.info(f'opening gencode annotations: {gencode}')
        cdef vector[NamedTx] transcripts
        cdef Gene curr
        if gencode is not None:
            transcripts = open_gencode(str(gencode).encode('utf8'), coding_only)
            for x in transcripts:
                symbol = x.symbol.decode('utf8')
                if symbol not in self.genes:
                    self.genes[symbol] = Gene(symbol.encode('utf8'))
                curr = self.genes[symbol]
                curr.add_tx(x.tx, x.is_principal)
                self.genes[symbol] = curr
        self._sort()
    
    def _sort(self):
        ''' index by starts and ends, to speed finding genes in a region
        '''
        for symbol in self.genes:
            gene = self.genes[symbol]
            chrom = gene.chrom.encode('utf8')
            symbol = symbol.encode('utf8')
            
            # ensure the chromosome is present
            if self.starts.count(chrom) == 0:
                self.starts[chrom] = []
            if self.ends.count(chrom) == 0:
                self.ends[chrom] = []
            
            self.starts[chrom].push_back(GenePoint(gene.start, symbol))
            self.ends[chrom].push_back(GenePoint(gene.end, symbol))
        
        # sort start and end coords by position
        for x, values in self.starts:
            self.starts[x] = sorted(values, key=lambda x: x['pos'])
        for x, values in self.ends:
            self.ends[x] = sorted(values, key=lambda x: x['pos'])
    
    def __repr__(self):
        return f'Gencode(n_genes={len(self)})'
    def __len__(self):
        return len(self.genes)
    def __getitem__(self, symbol):
        return self.genes[symbol]
    def __iter__(self):
        for x in self.genes:
            yield x
    
    def add_gene(self, gene):
        ''' add another gene to the Gencode object
        '''
        if gene.symbol not in self.genes:
            self.genes[gene.symbol] = gene
        self._sort()
    
    def nearest(self, str chrom, int pos):
        ''' find the nearest gene to a genomic chrom, pos coordinate
        '''
        chrom = f'chr{chrom}' if not chrom.startswith('chr') else chrom
        _chrom = chrom.encode('utf8')
        
        # first, account for any overlapping genespython
        overlaps = self.in_region(chrom, pos-1, pos+1)  # NOTE: possibly fix?
        if len(overlaps) > 0:
            # if we have > 0 prioritise if the position is in the CDS
            cds_overlaps = [x for x in overlaps if x.in_any_tx_cds(pos)]
            if len(cds_overlaps) > 0:
                overlaps = cds_overlaps
            # prioritise the gene with longest CDS (in the canonical tx)
            txs = [x.canonical for x in overlaps]
            lengths = [x.get_coding_distance(x.get_cds_end())['pos'] for x in txs]
            idx = lengths.index(max(lengths))
            return overlaps[idx]
        
        # no overlaps observed, look for the nearest upstream or downstream gene
        cdef GenePoint site = GenePoint(pos, b'A');
        cdef int i = lower_bound(self.starts[_chrom].begin(), self.starts[_chrom].end(), site, CompFunc) - self.starts[_chrom].begin()
        cdef int j = upper_bound(self.ends[_chrom].begin(), self.ends[_chrom].begin(), site, CompFunc) - self.ends[_chrom].begin()
        
        i = min(i, self.starts[_chrom].size() - 1)
        j = min(j, self.starts[_chrom].size() - 1)
        
        downstream = self[self.starts[_chrom][i].symbol.decode('utf8')]
        upstream = self[self.starts[_chrom][j].symbol.decode('utf8')]
        
        if upstream.distance(chrom, pos) <= downstream.distance(chrom, pos):
            return upstream
        else:
            return downstream
    
    def in_region(self, str _chrom, int start, int end, int max_window=2500000):
        ''' find genes within a genomic region
        
        Args:
            chrom: chromosome to search on
            start: start position of region
            end: end position of region
            max_window: some genes encapsulate the region, which means we have 
                to account for gene lengths of up to 2.3 Mb in the human genome.
                This permits extra search space in other organisms.
        
        Returns:
            list of Gene objects
        '''
        symbols = _in_region(_chrom.encode('utf8'), start, end, self.starts,
            self.ends, max_window)
        return [self[x.decode('utf8')] for x in symbols]
    
    def __exit__(self):
        ''' cleanup at exit
        '''
        global __genome_
        __genome_.close()
        __genome_ = None
