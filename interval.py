""" class to hold interval objects (such as genes, transcript regions etc)
"""

from sequence_methods import SequenceMethods

class Interval(SequenceMethods):
    """ class to hold intervals for gene regions, including CDS positions
    
    Uses BED gene and exon coordinates to define positions, and check whether
    submitted positions lie within the interval, or inside the CDS regions. 
    Also can determine distance between two positions in the CDS region.
    """
    
    # def __init__(self, bed_line):
    def __init__(self, transcript_id, start, end, strand, chrom, exon_ranges, cds_ranges):
        """ initiate the class with a line from a BED file
        
        Args:
            bed_line: line from bed file (11 column format)
        """
        
        self.chrom = chrom
        self.start = long(start)
        self.end = long(end)
        self.name = transcript_id
        self.strand = strand
        self.cds_start = ""
        self.cds_end = ""
        
        self.exons = self.__add_exons__(exon_ranges)
        self.cds = self.__add_cds__(cds_ranges)
        
    def __add_exons__(self, exon_ranges):
        """ add exon positions from the BED line
        """
        
        exon_ranges = sorted(exon_ranges)
        
        corrected_exon_ranges = []
        for start, end in exon_ranges:
            corrected_exon_ranges.append((start, end + 1))
        
        return corrected_exon_ranges
    
    def __add_cds__(self, cds_ranges):
        """ add cds positions from the exon positions and the CDS start and end.
        """
        
        cds = sorted(cds_ranges)
        fixed_cds = []
        for start, end in cds:
            fixed_cds.append((start, end + 1))
        
        cds = fixed_cds
        
        self.cds_start = long(cds[0][0])
        self.cds_end = long(cds[-1][-1])
        
        # sometimes the cds start and end lie outside the exons. We find the 
        # offset to the closest exon boundary, and use that to place the 
        # position within the up or downstream exon
        if not self.in_exons(self.cds_start):
            self.cds_start = self.fix_out_of_exon_cds_boundary(self.cds_start)
        if not self.in_exons(self.cds_end):
            self.cds_end = self.fix_out_of_exon_cds_boundary(self.cds_end)
        
        return cds
        
    def fix_out_of_exon_cds_boundary(self, position):
        """ sometimes the cds start and end lie outside the exons. We adjust the
        position to the prior or following exon
        """
        
        (start, end) = self.find_closest_exon(position)
        
        start_dist = abs(position - start)
        end_dist = abs(position - end)
        
        # print "CDS out of bounds:", self.get_name(), "by", min(start_dist, end_dist)
        
        # if we are closer to the start, then we go back an exon
        if start_dist < end_dist:
            before = self.get_exon_containing_position(start, self.exons) - 1
            position = self.exons[before][1] - start_dist 
        # if we are closer to the end, then we go forward an exon
        elif start_dist > end_dist:
            after = self.get_exon_containing_position(end, self.exons) + 1
            position = self.exons[after][0] + end_dist 
        
        return position
        
    def convert_chrom_to_int(self, chrom):
        """ converts a chromosome string to an int (if possible) for sorting.
        
        Args: 
            chrom: string (eg "1", "2", "3" ... "22", "X", "Y")
        
        Returns:
            int value of chrom
        """
        
        # set the integer values of the sex chromosomes
        chrom_dict = {"X": 23, "CHRX": 23, "Y": 24, "CHRY": 24, "MT": 25, \
            "CHRMT": 25}
        
        try:
            chrom = int(chrom)
        except ValueError:
            chrom = chrom_dict[chrom.upper()]
        
        return chrom
    
    def get_chrom(self):
        """ returns the chrom name for the BED interval
        """
        return self.chrom
    
    def get_name(self):
        """ returns the transcript ID (aka BED name) for the BED interval
        """
        return self.name
    
    def get_start(self):
        """ returns the start position for the interval (note not the CDS start)
        """
        
        return self.start
    
    def get_end(self):
        """ returns the end position for the interval (note not the CDS end)
        """
        
        return self.end
    
    def get_cds_start(self):
        """ returns the cds start position for the interval
        """
        
        return self.cds_start
    
    def get_cds_end(self):
        """ returns the cds end position for the interval
        """
        
        return self.cds_end
    
    def __eq__(self, other):
        return self.get_chrom() == other.get_chrom() and \
            self.get_start() == other.get_start() and \
            self.get_end() == other.get_end()
            
    def __ne__(self, other):
        return self.get_chrom() != other.get_chrom() or \
            self.get_start() != other.get_start() or \
            self.get_end() != other.get_end()
    
    def __lt__(self, other):
        return (self.convert_chrom_to_int(self.get_chrom()), self.get_start(), \
                self.get_end()) < \
                (other.convert_chrom_to_int(other.get_chrom()), \
                other.get_start(), other.get_end())
    
    def __repr__(self):
        return "{0}({1}:{2}-{3})".format(self.__class__.__name__, \
            self.get_chrom(), self.get_start(), self.get_end())
    
    def __str__(self):
        return "\t".join(self.line)
    
    def __hash__(self):
        return hash((self.get_chrom(), self.get_start(), self.get_end()))
    
    def contains_position(self, position):
        """ check if a position lies within the interval boundaries
        
        Args:
            position: nucleotide position in bp
        
        Returns:
            True/False for whether the position is within the gene
        """
        
        return self.start <= long(position) <= self.end
    
    def in_exons(self, position):
        """ determines if a nucleotide position lies within the exons
        """
        
        position = long(position)
        
        for start, end in self.exons:
            if start <= position <= end:
                return True
        
        return False
    
    def find_closest_exon(self, position):
        """ finds the positions of the exon closest to a position
        """
        
        max_distance = 999999999
        ref_start = 0
        ref_end = 0
        
        for start, end in self.exons:
            start_dist = abs(start - position)
            end_dist = abs(end - position)
            
            if start_dist < max_distance:
                ref_start = start
                ref_end = end
                max_distance = start_dist
            
            if end_dist < max_distance:
                ref_start = start
                ref_end = end
                max_distance = end_dist
        
        return (ref_start, ref_end)
    
    def in_coding_region(self, position):
        """ determines if a nucleotide position lies within the coding region
        
        Args:
            position:
        
        Returns:
            True/False for whether the position is within the coding region
        """
        
        position = long(position)
        
        for start, end in self.cds:
            if start <= position <= end:
                return True
        
        return False
    
    def get_exon_containing_position(self, position, exon_ranges):
        """ find the exon number for a position within an exon
        
        Args:
            position: chromosomal nucleotide position
        
        Returns:
            number of exon containing the position
        
        Raises:
            AssertionError if position not in coding region
        """
        
        position = long(position)
        
        exon_num = 0
        for start, end in exon_ranges:
            if start <= position <= end:
                return exon_num
            exon_num += 1
        
        return None
    
    def get_coding_distance(self, pos_1, pos_2):
        """ find the distance in coding bp between two chr positions
        
        Args:
            pos_1: first chromosome position
            pos_2: second chromosome position
        
        Returns:
            distance in coding sequence base pairs between two positions
        """
        
        
        # make sure that the positions are within the coding region, otherwise
        # there's no point trying to calculate the coding distance
        assert self.in_coding_region(pos_1)
        assert self.in_coding_region(pos_2)
        
        min_pos = min(pos_1, pos_2)
        max_pos = max(pos_1, pos_2)
        
        exon_1 = self.get_exon_containing_position(min_pos, self.cds)
        exon_2 = self.get_exon_containing_position(max_pos, self.cds)
        
        if exon_1 == exon_2:
            distance = max_pos - min_pos
        else:
            exon_1_distance = self.cds[exon_1][1] - min_pos
            exon_2_distance = max_pos - self.cds[exon_2][0]
            
            distance = exon_1_distance + exon_2_distance
            for exon in range(exon_1 + 1, exon_2):
                distance += self.cds[exon][1] - self.cds[exon][0]
        
        return distance
    
        
        
