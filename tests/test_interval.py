""" class to test the Interval class
"""

import unittest

from denovonear.interval import Interval

class TestIntervalPy(unittest.TestCase):
    """ unit test the Interval class
    """
    
    def setUp(self):
        """ construct an interval object for unit tests
        """
        
        chrom = 1
        start = 1000
        end = 2000
        name = "TEST"
        strand = "+"
        exons = [(1000, 1200), (1800, 2000)]
        cds = [(1100, 1200), (1800, 1900)]
        
        self.gene = Interval(name, start, end, strand, chrom, exons, cds)
    
    def test___add_exons__(self):
        """ test that __add_exons__() works correctly
        """
        
        exons = [(0, 200), (800, 1000)]
        cds = [(1100, 1200), (1800, 1900)]
        self.assertEqual(self.gene.__add_exons__(exons, cds), [(0, 200), (800, 1000)])
        
        self.gene.strand = "-"
        self.assertEqual(self.gene.__add_exons__(exons, cds), [(-1, 199), (799, 999)])
        
    
    def test___add_exons___missing_exon(self):
        """ test that __add_exons__() works correctly when we lack coordinates
        """
        
        # also test that we can determine the exons when we don't have any given
        # for a transcript with a single CDS region
        self.gene.strand = "+"
        exons = []
        cds = [(1100, 1200)]
        self.assertEqual(self.gene.__add_exons__(exons, cds), [(1000, 2000)])
        
        # check that missing exons, but 2+ CDS regions raises an error.
        cds = [(1100, 1200), (1300, 1400)]
        with self.assertRaises(ValueError):
            self.gene.__add_exons__(exons, cds)
        
        # if the CDS start matches the gene start, check that the gene is
        # extended by 1 bp
        self.gene.start = 1000
        self.gene.end = 2000
        cds = [(1000, 1200)]
        self.assertEqual(self.gene.__add_exons__(exons, cds), [(999, 2000)])
        
        # if the CDS end matches the gene end, check that the gene is 
        # extended by 1 bp
        self.gene.start = 1000
        self.gene.end = 2000
        cds = [(1100, 2000)]
        self.assertEqual(self.gene.__add_exons__(exons, cds), [(1000, 2001)])
    
    def test___add_cds__(self):
        """ test that __add_cds__() works correctly
        """
        
        self.gene.exons = [(0, 200), (800, 1000)]
        
        cds = [(100, 200), (800, 900)]
        
        # check that positions on the minus strand are adjusted back one base
        self.gene.strand = "-"
        self.assertEqual(self.gene.__add_cds__(cds), [(99, 199), (799, 899)])
        
        self.gene.strand = "+"
        self.assertEqual(self.gene.__add_cds__(cds), cds)
        
        # check that CDS ends outside an exon are corrected
        self.gene.exons = [(0, 200), (300, 400), (800, 1000)]
        cds = [(100, 200), (300, 402)]
        self.assertEqual(self.gene.__add_cds__(cds), [(100, 200), (300, 400), (800, 802)])
        
        cds = [(298, 400), (800, 1000)]
        self.assertEqual(self.gene.__add_cds__(cds), [(198, 200), (300, 400), (800, 1000)])
    
    def test_fix_out_of_exon_cds_boundary(self):
        """ test that fix_out_of_exon_cds_boundary() works correctly
        """
        
        self.gene.exons = [(1100, 1200), (1300, 1400), (1800, 1900)]
        
        self.assertEqual(self.gene.fix_out_of_exon_cds_boundary(1295), (1195, 1200))
        self.assertEqual(self.gene.fix_out_of_exon_cds_boundary(1205), (1300, 1305))
        
        self.assertEqual(self.gene.fix_out_of_exon_cds_boundary(1402), (1800, 1802))
        self.assertEqual(self.gene.fix_out_of_exon_cds_boundary(1798), (1398, 1400))
        
        # raise an error if the position is within the exons
        with self.assertRaises(AssertionError):
            self.gene.fix_out_of_exon_cds_boundary(1105)
        
    
    def test_convert_chrom_to_int(self):
        """ test that convert_chrom_to_int() works correctly
        """
        
        self.assertEqual(self.gene.convert_chrom_to_int("1"), 1)
        self.assertEqual(self.gene.convert_chrom_to_int("2"), 2)
        self.assertEqual(self.gene.convert_chrom_to_int("X"), 23)
        self.assertEqual(self.gene.convert_chrom_to_int("chrx"), 23)
        self.assertEqual(self.gene.convert_chrom_to_int("chrX"), 23)
        self.assertEqual(self.gene.convert_chrom_to_int("y"), 24)
        self.assertEqual(self.gene.convert_chrom_to_int("chry"), 24)
        self.assertEqual(self.gene.convert_chrom_to_int("mt"), 25)
        self.assertEqual(self.gene.convert_chrom_to_int("chrmt"), 25)
        self.assertRaises(KeyError, self.gene.convert_chrom_to_int, "Z")
    
    def test___add__(self):
        """ test that __add__() works correctly
        """
        
        exons = [(10, 20), (50, 60), (90, 100)]
        a = Interval("a", 10, 100, "+", "1", exons, [(55, 60), (90, 100)])
        b = Interval("b", 10, 100, "+", "1", exons, [(50, 60), (90, 95)])
        c = Interval("c", 10, 100, "+", "1", [(45, 65)], [(45, 65)])
        d = Interval("d", 10, 100, "+", "1", [(30, 40)], [(30, 40)])
        
        # check that adding two Intervals gives the union of CDS regions
        self.assertEqual((a + b).cds, [(50, 60), (90, 100)])
        self.assertEqual((a + c).cds, [(45, 65), (90, 100)])
        
        # check that addition is reversible
        self.assertEqual((c + a).cds, [(45, 65), (90, 100)])
        
        # check that adding previously unknown exons works
        self.assertEqual((a + d).cds, [(30, 40), (55, 60), (90, 100)])
    
    def test_region_overlaps_cds(self):
        """ check that region_overlaps_cds() works correctly
        """
        
        # the cds regions are at [(1100, 1200), (1800, 1900)], so check regions
        # that do and do not intersect with those
        
        self.assertTrue(self.gene.region_overlaps_cds((1050, 1150)))
        
        # check exons surrounding, and within the genes exons
        self.assertTrue(self.gene.region_overlaps_cds((1050, 1250)))
        self.assertTrue(self.gene.region_overlaps_cds((1150, 1160)))
        
        # check that non overlapping region fails
        self.assertFalse(self.gene.region_overlaps_cds((1050, 1090)))
        
        # check the boundaries of the exons
        self.assertTrue(self.gene.region_overlaps_cds((1050, 1100)))
        self.assertFalse(self.gene.region_overlaps_cds((1050, 1099)))
    
    def test_in_exons(self):
        """ test that in_exons() works correctly
        """
        
        self.gene.exons = [(1000, 1200), (1800, 2000)]
        
        # check for positions inside the exon ranges
        self.assertTrue(self.gene.in_exons(1000))
        self.assertTrue(self.gene.in_exons(1001))
        self.assertTrue(self.gene.in_exons(1200))
        self.assertTrue(self.gene.in_exons(1800))
        self.assertTrue(self.gene.in_exons(1801))
        self.assertTrue(self.gene.in_exons(1999))
        self.assertTrue(self.gene.in_exons(2000))
        
        # check positions outside the exon ranges
        self.assertFalse(self.gene.in_exons(999))
        self.assertFalse(self.gene.in_exons(1201))
        self.assertFalse(self.gene.in_exons(1799))
        self.assertFalse(self.gene.in_exons(2001))
        self.assertFalse(self.gene.in_exons(-1100))
    
    def test_find_closest_exon(self):
        """ test that find_closest_exon() works correctly
        """
        
        exon_1 = (1000, 1200)
        exon_2 = (1800, 2000)
        self.gene.exons = [exon_1, exon_2]
        
        # find for positions closer to the first exon
        self.assertEqual(self.gene.find_closest_exon(0), exon_1)
        self.assertEqual(self.gene.find_closest_exon(999), exon_1)
        self.assertEqual(self.gene.find_closest_exon(1000), exon_1)
        self.assertEqual(self.gene.find_closest_exon(1100), exon_1)
        self.assertEqual(self.gene.find_closest_exon(1200), exon_1)
        self.assertEqual(self.gene.find_closest_exon(1201), exon_1)
        self.assertEqual(self.gene.find_closest_exon(1500), exon_1)
        
        # find for positions closer to the second exon
        self.assertEqual(self.gene.find_closest_exon(1501), exon_2)
        self.assertEqual(self.gene.find_closest_exon(1799), exon_2)
        self.assertEqual(self.gene.find_closest_exon(1800), exon_2)
        self.assertEqual(self.gene.find_closest_exon(1900), exon_2)
        self.assertEqual(self.gene.find_closest_exon(2000), exon_2)
        self.assertEqual(self.gene.find_closest_exon(2001), exon_2)
        self.assertEqual(self.gene.find_closest_exon(10000), exon_2)
    
    def test_in_coding_region(self):
        """ test that in_coding_region() works correctly
        """
        
        self.gene.cds = [(1100, 1200), (1800, 1900)]
        
        # check for positions inside the exon ranges
        self.assertTrue(self.gene.in_coding_region(1100))
        self.assertTrue(self.gene.in_coding_region(1101))
        self.assertTrue(self.gene.in_coding_region(1200))
        self.assertTrue(self.gene.in_coding_region(1800))
        self.assertTrue(self.gene.in_coding_region(1801))
        self.assertTrue(self.gene.in_coding_region(1899))
        self.assertTrue(self.gene.in_coding_region(1900))
        
        # check positions outside the exon ranges
        self.assertFalse(self.gene.in_coding_region(1099))
        self.assertFalse(self.gene.in_coding_region(1201))
        self.assertFalse(self.gene.in_coding_region(1799))
        self.assertFalse(self.gene.in_coding_region(1901))
        self.assertFalse(self.gene.in_coding_region(-1100))
    
    def test_get_exon_containing_position(self):
        """ test that get_exon_containing_position() works correctly
        """
        
        
        
        exons = [(1000, 1200), (1800, 2000)]
        
        self.assertEqual(self.gene.get_exon_containing_position(1000, exons), 0)
        self.assertEqual(self.gene.get_exon_containing_position(1200, exons), 0)
        self.assertEqual(self.gene.get_exon_containing_position(1800, exons), 1)
        self.assertEqual(self.gene.get_exon_containing_position(2000, exons), 1)
        
        # raise an error if the position isn't within the exons
        with self.assertRaises(ValueError):
            self.gene.get_exon_containing_position(2100, exons)
    
    def test_get_coding_distance(self):
        """ test that get_coding_distance() works correctly
        """
        
        self.gene.cds = [(1100, 1200), (1800, 1900)]
        
        # raise an error for positions outside the CDS
        with self.assertRaises(AssertionError):
            self.gene.get_coding_distance(1000, 1100)
        with self.assertRaises(AssertionError):
            self.gene.get_coding_distance(1100, 1300)
        
        # zero distance between a site and itself
        self.assertEqual(self.gene.get_coding_distance(1100, 1100), 0)
        
        # within a single exon, the distance is between the start and end
        self.assertEqual(self.gene.get_coding_distance(1100, 1200), 100)
        
        # if we traverse exons, the distance bumps up at exon boundaries
        self.assertEqual(self.gene.get_coding_distance(1100, 1800), 101)
        
        # check full distance across gene
        self.assertEqual(self.gene.get_coding_distance(1100, 1900), 201)
        
        # check that the distance bumps up for each exon boundary crossed
        self.gene.cds = [(1100, 1200), (1300, 1400), (1800, 1900)]
        self.assertEqual(self.gene.get_coding_distance(1100, 1900), 302)
    
    def test_convert_chr_pos_to_cds_positions(self):
        """ test that convert_chr_pos_to_cds_positions() works correctly
        """
        
        self.gene.cds = [(1100, 1200), (1800, 1900)]
        self.gene.strand = "+"
        
        # note that all of these chr positions are 0-based (ie pos - 1)
        self.assertEqual(self.gene.convert_chr_pos_to_cds_positions(1099), 0)
        self.assertEqual(self.gene.convert_chr_pos_to_cds_positions(1100), 1)
        self.assertEqual(self.gene.convert_chr_pos_to_cds_positions(1199), 100)
        
        # check that outside exon boundaries gets the closest exon position, if
        # the variant is close enough
        self.assertEqual(self.gene.convert_chr_pos_to_cds_positions(1200), 100)
        self.assertEqual(self.gene.convert_chr_pos_to_cds_positions(1201), 100)
        self.assertEqual(self.gene.convert_chr_pos_to_cds_positions(1798), 101)
        self.assertEqual(self.gene.convert_chr_pos_to_cds_positions(1799), 101)
        
        # check that sites sufficiently distant from an exon raise an error, or
        # sites upstream of a gene, just outside the CDS, but within an exon
        with self.assertRaises(ValueError):
            self.gene.convert_chr_pos_to_cds_positions(1215)
        with self.assertRaises(ValueError):
            self.gene.convert_chr_pos_to_cds_positions(1098)
        
        # check that sites in a different exon are counted correctly
        self.assertEqual(self.gene.convert_chr_pos_to_cds_positions(1799), 101)
        
        # check that sites on the reverse strand still give the correct CDS
        self.gene.strand = "-"
        self.assertEqual(self.gene.convert_chr_pos_to_cds_positions(1899), 0)
        self.assertEqual(self.gene.convert_chr_pos_to_cds_positions(1890), 9)
        self.assertEqual(self.gene.convert_chr_pos_to_cds_positions(1799), 100)
        self.assertEqual(self.gene.convert_chr_pos_to_cds_positions(1790), 100)
        
        self.assertEqual(self.gene.convert_chr_pos_to_cds_positions(1200), 101)
        
        
    
    


# if __name__ == '__main__':
#     unittest.main()
