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

cdef class Transcript:
    def __cinit__(self, transcript_id, chrom, start, end, strand):
        transcript_id = transcript_id.encode('utf8')
        chrom = chrom.encode('utf8')
        self.thisptr = new Tx(transcript_id, chrom, start, end, ord(strand))
    def __dealloc__(self):
        del self.thisptr
    
    def __repr__(self):
        return "{0}({1} {2}:{3}-{4})".format(self.__class__.__name__,
            self.thisptr.get_name(), self.thisptr.get_chrom(),
            self.thisptr.get_start(), self.thisptr.get_end())
    
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
        
        self.thisptr.set_exons(exon_ranges, cds_ranges)
    
    
    def get_overlap(self, exon, regions):
        ''' find the union of ranges for overlapping regions
        
        Args:
            exon: dictionary of start and end positions for a single exon
            regions: list of start and end position dictionaries.
        
        Returns:
            list of start and end positions, or None if no matches found
        '''
        
        for region in regions:
            if (exon['start'] <= region['end'] and exon['end'] >= region['start']):
                return {'start': min(exon['start'], region['start']),
                    'end': max(exon['end'], region['end'])}
        
        return None
    
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
        for x in first:
            overlap = self.get_overlap(x, second)
            
            if overlap is not None:
                coords.append(overlap)
            else:
                coords.append(x)
        
        # include the second's regions that didn't overlap the first's regions
        for x in second:
            overlap = self.get_overlap(x, coords)
            
            if overlap is None:
                coords.append(x)
        
        coords = [ ( x['start'], x['end'] ) for x in coords ]
        
        return sorted(coords)
    
    def merge_genomic_seq(self, other):
        ''' merge the genomic sequence from two transcripts
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
        
        final = lag_gdna[not_lag.get_end() - lag.get_start() + lead_offset:]
        
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
        
        altered = Transcript('{}:{}'.format(self.get_name(), other.get_name()),
            self.get_chrom(), min(self.get_start(), other.get_start()),
            min(self.get_end(), other.get_end()), self.get_strand())
        
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
        
    def chrom_pos_to_cds(self, pos_1):
        return self.thisptr.chrom_pos_to_cds(pos_1)
    
    def get_position_on_chrom(self, pos):
        return self.thisptr.get_position_on_chrom(pos)
    
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
