""" class to test the Interval class
"""

import unittest
from interval import Interval

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
        self.assertEqual(self.gene.__add_exons__(exons), [(0, 200), (800, 1000)])
        
        self.gene.strand = "-"
        self.assertEqual(self.gene.__add_exons__(exons), [(-1, 199), (799, 999)])
    
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
    
    def test_contains_position(self):
        """ test that contains_position() works correctly
        """
        
        self.gene.start = 1000
        self.gene.end = 1100
        
        self.assertTrue(self.gene.contains_position(1000))
        self.assertTrue(self.gene.contains_position(1001))
        self.assertTrue(self.gene.contains_position(1100))
        self.assertFalse(self.gene.contains_position(999))
        self.assertFalse(self.gene.contains_position(1101))
        
        # also test that it works at chrom positions well beyond what we 
        # expect to encounter. This will be the only check that it can handle
        # large positions
        self.gene.start = 1000001000
        self.gene.end = 1000001100
        
        self.assertTrue(self.gene.contains_position(1000001000))
        self.assertTrue(self.gene.contains_position(1000001100))
        self.assertFalse(self.gene.contains_position(1000000999))
        self.assertFalse(self.gene.contains_position(1000001101))
    
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
    
    


if __name__ == '__main__':
    unittest.main()

