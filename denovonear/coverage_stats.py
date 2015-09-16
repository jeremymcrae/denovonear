""" allows for coverage-based adjustments to mutation rate probabilities
"""

import math
import pysam

class CoverageStats(object):
    """ class to obtain coverage information for a genomic region, and adjust 
    mutation rate probabilities given coverage information.
    """
    
    tabix = None
    chrom = None
    columns = {1: 4, 5: 5, 10: 6, 15: 7, 20: 8, 25: 9, 30: 10, 50: 11, 100: 12}
    coverage_dir = "ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.1/coverage"
    region_coverages = {}
    
    def set_coverage_dir(self, coverage_dir):
        """ set the location of the coverage files
        
        Args:
            coverage_dir: folder containing coverage file, can be ftp folder, or
                local directory
        """
        
        self.coverage_dir = coverage_dir
    
    def sanitise_chrom(self, chrom):
        """ make sure the chromosome string is appropriate
        
        Args:
            chrom: integer, or string for chromosome
        
        Returns:
            chromosome as string, without any "chr" prefix
        """
        
        return str(chrom).upper().lstrip("CHR")
    
    def set_chrom(self, chrom):
        """ set the chromosome to use
        
        There are a bunch of chromosome functions, since we have to connect to a
        tabix file for a single chromosome, but we will also have to switch
        between different chromosome files. Thus we try to handle this correctly
        
        Args:
            chrom: integer, or string for chromosome
        """
        
        self.chrom = self.sanitise_chrom(chrom)
    
    def get_chrom(self):
        return self.chrom
    
    def set_tabix(self):
        """ connect to the tabix indexed coverage dataset
        """
        
        # close the previous tabix file, if it is currently open
        if self.tabix is not None:
            self.tabix.close()
        
        path = "{0}/Panel.chr{1}.coverage.txt.gz".format(self.coverage_dir, self.get_chrom())
        
        try:
            self.tabix = pysam.Tabixfile(path)
        except IOError:
            self.tabix = None
    
    def set_coverages(self, chrom, start_pos, end_pos, min_coverage=10):
        """ find the proportion of samples exceeding a coverage in a region.
        
        The lines in the coverage files provided by ExAC consist of chromosome
        coordinates, along with several columns providing the proportion of
        samples with coverage exceeding defined coverage levels. We extract the
        proportion at a given coverage level for the bases in the region.
        
        Example of coverage file columns:
        
        #chrom  pos    mean  median  1       5     10    15    20    25    30    50    100
        1       12546  0.00  0.00    0.0002  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
        
        Creates a dictionary of proportion of samples with coverage exceeding 
        min_coverage, indexed by chromosome position. Note that this only
        extracts the chromosome positions that exist in the file, some sites
        do not have coverage data, so they will not be in the dictionary.
        
        Args:
            chrom: integer, or string for chromosome
            start_pos: initial position for the 
            end_pos: final position in the gene range
            min_coverage: the coverage level that we wish to check
        """
        
        self.region_coverages = {}
        
        assert min_coverage in self.columns
        assert type(start_pos) == int
        assert type(end_pos) == int
        assert end_pos >= start_pos
        
        # make sure we are grabbing data from the correct tabix file
        if self.sanitise_chrom(chrom) != self.get_chrom():
            self.set_chrom(chrom)
            self.set_tabix()
        
        # get an iterator to the correct lines
        try:
            rows = self.tabix.fetch(chrom, start_pos, end_pos)
        except AttributeError:
            rows = []
        
        # extract the coverage at each site within the lines
        for line in rows:
            line = line.strip().split("\t")
            position = int(line[1])
            coverage = float(line[self.columns[min_coverage]])
            self.region_coverages[position] = coverage
    
    def adjust_rate_for_coverage(self, position, rate):
        """ adjusts a mutation rate for the coverage level at a position
        
        Adjust a mutation rate in a similar manner to that described in the
        supplementary note of Nature Genetics 944-950 (2014) doi:10.1038/ng.3050
        They describe their method as:
        
        "For each base, we looked up the number of trios in which all members 
        had 10x coverage or greater and used that number to determine the 
        appropriate discount. For bases with almost all trios having 10x 
        coverage, the probability of mutation was not adjusted. However, as the
        number of trios with 10x coverage dropped, the probability of mutation 
        was multiplied by an adjustment factor in between 0.9 and 1."
        
        I'm instead using the proportion of total samples (rather than trios)
        with coverage greater than 10X coverage, and I'm defaulting to the ExAC
        coverage datasets rather than determining the project specific figure.
        
        Args:
            position: integer for chromosome nucleotide position
            rate: mutation rate for the nucleotide position
        
        Returns:
            
        """
        
        if position in self.region_coverages:
            adjustment_factor = 0.9 + (0.1 * self.region_coverages[position])
            rate *= adjustment_factor
        
        return rate
        
        
        