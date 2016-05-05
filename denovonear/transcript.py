""" class to hold transcript objects (such as genes, transcript regions etc)
"""

import copy

from denovonear.transcript_sequence_methods import SequenceMethods

class Transcript(SequenceMethods):
    """ Class to define transcript regions.
    
    Uses gene and exon coordinates to define positions, and check whether
    submitted positions lie within the transcript, or inside the coding sequence.
    Also can determine distance between two positions in coding distance between
    two positions in the CDS.
    """
    
    def __init__(self, transcript_id, start, end, strand, chrom, exon_ranges, cds_ranges):
        """ initiate the class with details for a transcript
        
        Args:
            transcript_id: ensembl transcript ID
            start: nucleotide position at 5' position (on + strand of chrom)
            end: nucleotide position at 3' position (on + strand of chrom)
            strand: strand that the transcript is on ("+" or "-")
            chrom: string to indicate the chromosome that the transcript is on
            exon_ranges: list of (start, end) tuples that contain all the exon
                positions for the transcript.
            cds_ranges: list of (start, end) tuples that contain all the coding
                sequence positions for the transcript.
        """
        
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = transcript_id
        self.strand = strand
        self.cds_min = ""
        self.cds_max = ""
        
        self.exons = self.__add_exons__(exon_ranges, cds_ranges)
        self.cds = self.__add_cds__(cds_ranges)
        
    def __add_exons__(self, exon_ranges, cds):
        """ add exon positions
        """
        
        exon_ranges = sorted(exon_ranges)
        
        # If the transcript lacks exon coordinates and only has a single CDS
        # region, then if the CDS region fits within the gene range, make a
        # single exon, using the transcript start and end. This prevents issues
        # when adding the CDS coordinates for the transcript.
        if exon_ranges == []:
            if len(cds) == 1 and min(cds[0]) >= self.start and max(cds[0]) <= self.end:
                exon_ranges = [(self.start, self.end)]
            else:
                # We can figure out the exon coordinates given a single CDS
                # region, but we can't do that given multiple CDS regions. I
                # haven't hit this yet, deal with it when it happens.
                raise ValueError("{} lacks exon coordinates".format(self.name))
        
        return exon_ranges
    
    def __add_cds__(self, cds_ranges):
        """ add cds positions from the exon positions and the CDS start and end.
        """
        
        cds = sorted(cds_ranges)
        
        self.cds_min = int(cds[0][0])
        self.cds_max = int(cds[-1][-1])
        
        # sometimes the cds start and end lie outside the exons. We find the
        # offset to the closest exon boundary, and use that to place the
        # position within the up or downstream exon, and correct the initial
        # cds boundaries
        if not self.in_exons(self.cds_min):
            (start, end) = self.fix_out_of_exon_cds_boundary(self.cds_min)
            cds[0] = (cds[0][0] + abs(end - start), cds[0][1])
            self.cds_min = start
            cds.insert(0, (start, end))
        if not self.in_exons(self.cds_max):
            (start, end) = self.fix_out_of_exon_cds_boundary(self.cds_max)
            cds[-1] = (cds[-1][0], cds[-1][1] - abs(end - start))
            self.cds_max = end
            cds.append((start, end))
        
        return cds
        
    def fix_out_of_exon_cds_boundary(self, position):
        """ fix cds start and ends which occur outside exons.
        
        The main cause is when the stop codon overlaps an exon boundary. In
        some of these cases, the CDS end is placed within the adjacent intron
        sequence (if that allows stop site sequence). To fix this, we adjust the
        position to the offset distance within the adjacent exon.
        
        Args:
            position: chromosome position
        
        Returns:
            adjusted chromosome position
        """
        
        assert not self.in_exons(position)
        
        (start, end) = self.find_closest_exon(position)
        
        start_dist = abs(position - start)
        end_dist = abs(position - end)
        
        # if we are closer to the start, then we go back an exon
        if start_dist < end_dist:
            before = self.get_exon_containing_position(start, self.exons) - 1
            end = self.exons[before][1]
            start = end - start_dist
        # if we are closer to the end, then we go forward an exon
        elif start_dist > end_dist:
            after = self.get_exon_containing_position(end, self.exons) + 1
            start = self.exons[after][0]
            end = start + end_dist
        
        return (start, end)
        
    def chrom_to_int(self, chrom):
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
        """ returns the chrom name for the transcript
        """
        return self.chrom
    
    def get_strand(self):
        return self.strand
    
    def get_name(self):
        """ returns the transcript ID for the transcript
        """
        return self.name
    
    def get_start(self):
        """ returns the start position for the transcript (note not the CDS start)
        """
        
        return self.start
    
    def get_end(self):
        """ returns the end position for the transcript (note not the CDS end)
        """
        
        return self.end
    
    def get_cds_start(self):
        """ returns the cds start position for the transcript
        """
        
        if self.get_strand() == "+":
            return self.cds_min
        elif self.get_strand() == "-":
            return self.cds_max
        else:
            raise ValueError("unknown strand type {0}".format(self.get_strand()))
    
    def get_cds_end(self):
        """ returns the cds end position for the transcript
        """
        
        if self.get_strand() == "+":
            return self.cds_max
        elif self.get_strand() == "-":
            return self.cds_min
        else:
            raise ValueError("unknown strand type {0}".format(self.get_strand()))
    
    def __eq__(self, other):
        return self.get_chrom() == other.get_chrom() and \
            self.get_start() == other.get_start() and \
            self.get_end() == other.get_end()
            
    def __ne__(self, other):
        return self.get_chrom() != other.get_chrom() or \
            self.get_start() != other.get_start() or \
            self.get_end() != other.get_end()
    
    def __lt__(self, other):
        return (self.chrom_to_int(self.get_chrom()), self.get_start(), \
                self.get_end()) < \
                (other.chrom_to_int(other.get_chrom()), \
                other.get_start(), other.get_end())
    
    def __repr__(self):
        return "{0}({1} {2}:{3}-{4})".format(self.__class__.__name__, \
            self.get_name(), self.get_chrom(), self.get_start(), self.get_end())
    
    def __str__(self):
        return self.__repr__()
    
    def __hash__(self):
        return hash((self.get_chrom(), self.get_start(), self.get_end()))
    
    def __add__(self, other):
        """ combine the coding sequences of two Transcript objects
        
        When we determine the sites for sampling, occasioally we want to
        use sites from multiple alternative transcripts. We determine the sites
        for each transcript in turn, but mask the sites that have been collected
        in the preceeding transcripts. In order to be able to mask all previous
        trabnscripts, we need to combine the coding sequence of the transcripts
        as we move through them. This function performs the union of coding
        sequence regions between different transcripts.
        
        Args:
            other: a transcript to be combined with the current object.
        
        Returns:
            an altered instance of the class, where the coding sequence regions
            are the union of the coding regions of two Transcript objects. This
            disrupts the ability to get meaningingful sequence from the object,
            so don't try to extract sequence from the returned object.
        """
        
        cds_min = min(other.get_cds_start(), other.get_cds_end())
        cds_max = max(other.get_cds_start(), other.get_cds_end())
        
        altered = copy.deepcopy(self)
        
        for (start, end) in other.exons:
            # check non-coding exons
            if not other.in_coding_region(start) and not other.in_coding_region(end):
                # ignore exons that don't contain any coding sequence regions
                if not other.region_overlaps_cds((start, end)):
                    continue
                else:
                    # this is an exon where the start and end positions are
                    # untranslated, but the exon still contains a coding region.
                    # This is unlikely, but we'll just swap the start and end
                    # to the cds start and end positions.
                    start = cds_min
                    end = cds_max
            
            # if the exon overlaps the start of the CDS, use the CDS start site
            # rather than the exon start
            if not other.in_coding_region(start) and other.in_coding_region(end):
                start = cds_min
                
            # if the exon overlaps the end of the CDS, use the CDS end site
            # rather than the exon end
            if other.in_coding_region(start) and not other.in_coding_region(end):
                end = cds_max
            
            # if the CDS region is already in the current object,we don't need
            # to extend the current objects coordinates
            if altered.in_coding_region(start) and altered.in_coding_region(end):
                continue
            
            # figure out if the cds to be added intersects with the CDS of the
            # alternative transcript
            if not altered.region_overlaps_cds((start, end)):
                altered.cds.append((start, end))
            else:
                # if the transcripts overlap, find the overlapping exon and
                # extend its coordinates
                closest = altered.find_closest_exon(start, exons=altered.cds)
                exon_num = altered.get_exon_containing_position(closest[0], altered.cds)
                
                new_exon = (min(closest[0], start), max(closest[1], end))
                altered.cds[exon_num] = new_exon
        
        # and tidy up the CDS start and end positions (even though these have
        # become less meaningful now that the CDS positions are from a union of
        # two transcripts).
        altered.cds = sorted(altered.cds)
        altered.cds_min = altered.cds[0][0]
        altered.cds_max = altered.cds[-1][1]
        
        return altered
    
    def region_overlaps_cds(self, exon, exon_ranges=None):
        """ checks if a region intersects any part of the CDS
        """
        
        if exon_ranges is None:
            exon_ranges = self.cds
        
        # figure out if the cds to be added intersects with the CDS of the
        # alternative transcript
        overlap = False
        for cds_exon in exon_ranges:
            if self.has_overlap(cds_exon, exon):
                overlap = True
        
        return overlap
    
    def has_overlap(self, range_1, range_2):
        """ checks if two ranges overlaps
        
        Args:
            range_1: tuple of (start, end) positions
            range_2: tuple of (start, end) positions
        """
        
        assert len(range_1) == 2
        assert len(range_2) == 2
        
        return range_1[0] <= range_2[1] and range_1[1] >= range_2[0]
    
    def in_exons(self, position):
        """ determines if a nucleotide position lies within the exons
        """
        
        position = int(position)
        
        for start, end in self.exons:
            if start <= position <= end:
                return True
        
        return False
    
    def find_closest_exon(self, position, exons=None):
        """ finds the positions of the exon closest to a position
        
        Can optionally specify a different set of exon ranges to check from,
        such as the CDS positions.
        """
        
        if exons is None:
            exons = self.exons
        
        max_distance = 1e10
        ref_start = 0
        ref_end = 0
        
        for start, end in exons:
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
            position: chromosomal nucleotide position
        
        Returns:
            True/False for whether the position is within the coding region
        """
        
        position = int(position)
        
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
        """
        
        position = int(position)
        
        exon_num = 0
        for start, end in exon_ranges:
            if start <= position <= end:
                return exon_num
            exon_num += 1
        
        raise ValueError("shouldn't get here")
    
    def get_coding_distance(self, pos_1, pos_2):
        """ find the distance in coding bp between two chr positions
        
        Args:
            pos_1: first chromosome position
            pos_2: second chromosome position
        
        Returns:
            distance in coding sequence base pairs between two positions
        
        Raises:
            AssertionError if position not in coding region
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
            exon_1_distance = (self.cds[exon_1][1] - min_pos)
            exon_2_distance = (max_pos - self.cds[exon_2][0]) + 1
            
            distance = exon_1_distance + exon_2_distance
            for exon in range(exon_1 + 1, exon_2):
                distance += (self.cds[exon][1] - self.cds[exon][0]) + 1
        
        return distance
    
    def chrom_pos_to_cds(self, pos):
        """ returns a chromosome position as distance from CDS ATG start
        """
        
        # need to convert the de novo event positions into CDS positions
        cds_start = self.get_cds_start()
        
        try:
            return self.get_coding_distance(cds_start, pos)
        except AssertionError:
            # catch the splice site functional mutations
            exon = self.find_closest_exon(pos)
            
            distances = [ abs(x - pos) for x in exon ]
            distance = min(distances)
            site = exon[distances.index(distance)]
            
            # catch the few variants that lie near an exon, but that exon isn't
            # part of the coding sequence
            if not self.in_coding_region(site):
                raise AssertionError("Not near coding exon: {0} in transcript"
                    " {1}".format(pos, self.get_name()))
            
            # ignore positions outside the exons that are too distant from a boundary
            if distance >= 9:
                raise AssertionError("distance to exon ({0}) > 8 bp for {1}"
                    " in transcript {2}".format(distance, pos, self.get_name()))
            
            return self.get_coding_distance(cds_start, site)
