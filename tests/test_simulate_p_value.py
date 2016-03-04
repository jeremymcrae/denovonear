""" unit test the WeightedChoice class"""

from __future__ import division

import os
import unittest

from denovonear.weights import WeightedChoice
from denovonear.transcript import Transcript
from denovonear.site_specific_rates import SiteRates
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.simulate import get_p_value

class TestGetPValuePy(unittest.TestCase):
    """ unit test the simulation of p-values
    """
        
    def test_get_p_value(self):
        """ check that get_p_value() works correctly
        """
        
        chrom = "1"
        name = "TEST"
        strand = "+"
        start = 0
        end = 70
        exons = [(5, 58)]
        cds = [(5, 58)]
        transcript = Transcript(name, start, end, strand, chrom, exons, cds)
        
        cds = "ATGTGGGCTCCACCAGCAGCAATCATGGGATGGGCCCACCAAGAAGGTGGGTAA"
        gdna = "GGGGGATGTGGGCTCCACCAGCAGCAATCATGGGATGGGCCCACCAAGAAGGTGGGTAACCAGGCCCC"
        transcript.add_cds_sequence(cds)
        transcript.add_genomic_sequence(gdna)
        
        # load the sequence contect mutation rates, then assess each site in the
        # CDS.
        path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data',
            'forSanger_1KG_mutation_rate_table.txt')
        mut_dict = load_mutation_rates(path)
        rates = SiteRates(transcript, mut_dict)
        
        cq = 'missense'
        de_novos = [5, 5]
        iterations = 10000
        
        (obs, p_value) = get_p_value(transcript, rates, iterations, cq, de_novos)
        
        self.assertTrue(p_value < 0.04)
        self.assertEqual(obs, '0.0')
        
        de_novos = [5, 58]
        
        (obs, p_value) = get_p_value(transcript, rates, iterations, cq, de_novos)
        
        self.assertTrue(p_value == 1.0)
        self.assertEqual(obs, '53.0')
        
