'''
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
'''

from itertools import combinations

cdef class Transcript:
    def __cinit__(self, name, chrom, start, end, strand, exons=None,
            cds=None, sequence=None, offset=0):
        ''' construct a Transcript object
        
        Args:
            name: ID of the transcript
            start: position in bp at 5' edge of transcript (on + strand)
            end: position in bp at 3' edge of transcript (on + strand)
            exons: list of tuples defining start and end positions of exons
            cds: list of tuples defining start and end positions of CDS regions
            sequence: DNA sequence of genome region of the transcript.
            offset: how many base pairs the DNA sequence extends outwards
        '''
        
        name = name.encode('utf8')
        chrom = chrom.encode('utf8')
        self.thisptr = new Tx(name, chrom, start, end, ord(strand))
        
        if exons is not None and cds is not None:
            self.set_exons(exons, cds)
            self.set_cds(cds)
            
            if sequence is not None:
                self.add_genomic_sequence(sequence, offset)
    
    def __dealloc__(self):
        del self.thisptr
    
    def __repr__(self):
        
        exons = [ (x['start'], x['end']) for x in self.get_exons() ]
        cds = [ (x['start'], x['end']) for x in self.get_cds() ]
        seq = self.get_genomic_sequence()
        
        if len(seq) > 40:
            seq = seq[:20] + '...[{} bp]...'.format(len(seq) - 40) + seq[-20:]
        
        if exons == []:
            exons = None
        
        if cds == []:
            cds = None
        
        if seq == '':
            seq = None
        else:
            seq = '"' + seq + '"'
        
        return 'Transcript(name="{}", chrom="{}", start={}, end={}, strand="{}", ' \
            'exons={}, cds={}, sequence={}, offset={})'.format(self.get_name(),
                self.get_chrom(), self.get_start(), self.get_end(),
                self.get_strand(), exons, cds, seq, self.get_genomic_offset())
    
    def __str__(self):
        return self.__repr__()
    
    def __hash__(self):
        return hash((self.thisptr.get_chrom(), self.thisptr.get_start(),
            self.thisptr.get_end()))
    
    def __richcmp__(self, other, op):
        
        if op == 2:
            return self.__hash__() == other.__hash__()
        else:
            err_msg = "op {0} isn't implemented yet".format(op)
            raise NotImplementedError(err_msg)
    
    def set_exons(self, exon_ranges, cds_ranges):
        ''' add exon ranges
        
        Args:
            exon_ranges: a CDS position of the selected base.
        '''
        
        if len(cds_ranges) == 0:
            raise ValueError('CDS coordinates were not supplied.')
        
        self.thisptr.set_exons(exon_ranges, cds_ranges)
    
    def get_overlaps(self, exon, regions):
        ''' find all regions which overlap a given region
        '''
        return [ i for i, x in enumerate(regions) if
            exon['start'] <= x['end'] and exon['end'] >= x['start'] ]
    
    def insert_region(self, coords, region):
        ''' include a region into a list of regions
        
        To include a region, we have to check which pre-existing regions the new
        region overlaps, so any overlaps can be merged into a single region.
        
        Args:
            coords: list of {start: X, end: Y} dictionaries
            region: dict of {'start': X, 'end': Y} positions
        '''
        indices = self.get_overlaps(region, coords)
        overlaps = [ coords[i] for i in indices ]
        start = min( x['start'] for x in overlaps + [region] )
        end = max( x['end'] for x in overlaps + [region] )
        
        for i in sorted(indices, reverse=True):
            del coords[i]
        
        return coords + [{'start': start, 'end': end}]
    
    def merge_coordinates(self, first, second):
        ''' merge two sets of coordinates, to get the union of regions
        
        This uses an inefficient approach, looping over and over, but we won't
        need to perform this often.
        
        Args:
            first: list of {'start': x, 'end': y} dictionaries for first transcript
            second: list of {'start': x, 'end': y} dictionaries for second transcript
        
        Returns:
            list of [start, end] lists, sorted by position.
        '''
        coords = []
        for a, b in combinations(first + second, 2):
            if a['start'] <= b['end'] and a['end'] >= b['start']:
                region = {'start': min(a['start'], b['start']),
                    'end': max(a['end'], b['end'])}
                a, b = region, region
            
            coords = self.insert_region(coords, a)
            coords = self.insert_region(coords, b)
        
        return [ (x['start'], x['end']) for x in sorted(coords, key=lambda x: x['start']) ]
    
    def merge_genomic_seq(self, other):
        ''' merge the genomic sequence from two transcripts
        
        We have two transcripts A, and B. We need to get the contiguous sequence
        from the start of the first transcript on the chromosome to the end of
        the second transcript. The transcripts may or may not overlap. There are
        three scenarios we need to account for:
        
        
        overlap without   A   =================
        enveloping                           ================= B
        
        
            overlap       A   ==============================
        and envelop               ==================  B
        
        
        no overlap        A  ===============
                                               ===============  B
        
        I've called the transcript whose sequence is first along the chromosome
        as 'lead', and the transcript whose sequence is last as 'lag', and the
        converse as 'not_lead', and 'not_lag'. Note that in the envelope case,
        the lead transcript can also be the lag transcript.
        '''
        
        # make sure that the surrounding sequence is the same length in both
        # transcripts.
        # TODO: this could be worked around, by figuring the minimum offset length,
        # TODO: then trimming the respective DNA offset sequences to that length.
        assert self.get_genomic_offset() == other.get_genomic_offset()
        
        # figure out which transcripts hold the leading and lagging sections
        lead, not_lead = self, other
        if self.get_start() > other.get_start():
            lead, not_lead = other, self
        
        lag, not_lag = self, other
        if self.get_end() < other.get_end():
            lag, not_lag = other, self
        
        lead_offset = lead.get_genomic_offset()
        lead_gdna = lead.get_genomic_sequence()
        initial = lead_gdna[:not_lead.get_start() - lead.get_start() + lead_offset]
        
        if self.get_start() <= other.get_end() and self.get_end() >= other.get_start():
            intersect_start = not_lead.get_start() - lead.get_start() + lead_offset
            intersect_end = not_lag.get_end() - lead.get_start() + lead_offset
            intersect = lead_gdna[intersect_start:intersect_end]
        else:
            intersect = 'N' * (lag.get_start() - lead.get_end() - lead_offset * 2)
        
        lag_offset = lag.get_genomic_offset()
        lag_gdna = lag.get_genomic_sequence()
        
        # some transcripts overlap, but some do not. We need to find the position
        # where the lagging transcript takes over, which is either at the end of
        # not lagging transcript, or the start of the lagging transcript,
        # whichever is higher
        pos = max(not_lag.get_end(), lag.get_start())
        final = lag_gdna[pos - lag.get_start() + lead_offset:]
        
        return initial + intersect + final
    
    def __add__(self, other):
        """ combine the coding sequences of two Transcript objects
        
        When we determine the sites for sampling, occasioally we want to
        use sites from multiple alternative transcripts. We determine the sites
        for each transcript in turn, but mask the sites that have been collected
        in the preceeding transcripts. In order to be able to mask all previous
        trabnscripts, we need to combine the coding sequence of the transcripts
        as we move through them. This function performs the union of coding
        sequence regions between different transcripts.
        
        We do this outside of the c++ class, so as to be able to set up a
        Transcript object correctly.
        
        Args:
            other: a transcript to be combined with the current object.
        
        Returns:
            an altered instance of the class, where the coding sequence regions
            are the union of the coding regions of two Transcript objects. This
            disrupts the ability to get meaningingful sequence from the object,
            so don't try to extract sequence from the returned object.
        """
        
        # if we try transcript + None or None + transcript, return the original
        # transcript, rather than raising an error.
        if other is None:
            return self
        if self is None:
            return other
        
        altered = Transcript('{}:{}'.format(self.get_name(), other.get_name()),
            self.get_chrom(), min(self.get_start(), other.get_start()),
            max(self.get_end(), other.get_end()), self.get_strand())
        
        exons = self.merge_coordinates(self.get_exons(), other.get_exons())
        cds = self.merge_coordinates(self.get_cds(), other.get_cds())
        
        altered.set_exons(exons, cds)
        altered.set_cds(cds)
        
        if self.get_genomic_sequence() != "":
            altered.add_genomic_sequence(self.merge_genomic_seq(other), self.get_genomic_offset())
        
        return altered
    
    def set_cds(self, cds_ranges):
        ''' set CDS ranges
        
        Args:
            cds_ranges: a CDS position of the selected base.
        '''
        
        self.thisptr.set_cds(cds_ranges)
    
    def get_genomic_offset(self):
        return self.thisptr.get_genomic_offset()
    def get_exons(self):
        return self.thisptr.get_exons()
    def get_cds(self):
        return self.thisptr.get_cds()
    def get_name(self):
        return self.thisptr.get_name().decode('utf8')
    def get_chrom(self):
        return self.thisptr.get_chrom().decode('utf8')
    def get_start(self):
        return self.thisptr.get_start()
    def get_end(self):
        return self.thisptr.get_end()
    def get_strand(self):
        return chr(self.thisptr.get_strand())
    def get_cds_start(self):
        return self.thisptr.get_cds_start()
    def get_cds_end(self):
        return self.thisptr.get_cds_end()
    def fix_cds_boundary(self, pos):
        return self.thisptr.fix_cds_boundary(pos)
    
    def in_exons(self, position):
        ''' check if a site lies within the exon ranges
        
        Args:
            position: an integer-based chromosome position.
        '''
        
        return self.thisptr.in_exons(position)
    
    def find_closest_exon(self, position, exons=None):
        ''' finds the positions of the exon closest to a position
        
        Can optionally specify a different set of exon ranges to check from,
        such as the CDS positions.
        '''
        
        cdef vector[Region] temp
        cdef Region pair
        
        if exons is None:
            return self.thisptr.find_closest_exon(position)
        else:
            for (start, end) in exons:
                pair.start = start
                pair.end = end
                temp.push_back(pair)
            return self.thisptr.find_closest_exon(position, temp)
    
    def in_coding_region(self, position):
        return self.thisptr.in_coding_region(position)
    
    def get_exon_containing_position(self, position, exons):
        '''
        '''
        
        cdef vector[Region] regions
        cdef Region pair
        
        for (start, end) in exons:
            pair.start = start
            pair.end = end
            regions.push_back(pair)
        
        return self.thisptr.get_exon_containing_position(position, regions)
    
    def get_coding_distance(self, pos_1, pos_2):
        return self.thisptr.get_coding_distance(pos_1, pos_2)
        
    def chrom_pos_to_cds(self, pos):
        coords = self.thisptr.chrom_pos_to_cds(pos)
        
        return {'pos': coords.position, 'offset': coords.offset}
    
    def get_position_on_chrom(self, pos, offset=0):
        return self.thisptr.get_position_on_chrom(pos, offset)
    
    def get_codon_number_for_cds_position(self, pos):
        return self.thisptr.get_codon_number_for_cds_position(pos)
    
    def get_position_within_codon(self, pos):
        return self.thisptr.get_position_within_codon(pos)
    
    def add_cds_sequence(self, text):
        self.thisptr.add_cds_sequence(text.encode('utf8'))
        
    def get_cds_sequence(self):
        return self.thisptr.get_cds_sequence().decode('utf8')
    
    def add_genomic_sequence(self, text, offset=0):
        self.thisptr.add_genomic_sequence(text.encode('utf8'), offset)
    
    def get_genomic_sequence(self):
        return self.thisptr.get_genomic_sequence().decode('utf8')
    
    def reverse_complement(self, text):
        return self.thisptr.reverse_complement(text).decode('utf8')
    
    def get_centered_sequence(self, pos, length=3):
        return self.thisptr.get_centered_sequence(pos, length).decode('utf8')
    
    def get_codon_sequence(self, pos):
        return self.thisptr.get_codon_sequence(pos).decode('utf8')
    
    def translate(self, text):
        return self.thisptr.translate(text.encode('utf8')).decode('utf8')
    
    def get_codon_info(self, pos):
        codon = dict(self.thisptr.get_codon_info(pos))
        
        if codon['codon_number'] == -1:
            codon['codon_number'] = None
            codon['intra_codon'] = None
            codon['codon_seq'] = None
            codon['initial_aa'] = None
        
        if codon['codon_seq'] is not None:
            codon['codon_seq'] = codon['codon_seq'].decode('utf8')
        
        if codon['initial_aa'] is not None:
            codon['initial_aa'] = codon['initial_aa'].decode('utf8')
        
        return codon
    
    def get_boundary_distance(self, pos):
        return self.thisptr.get_boundary_distance(pos)
