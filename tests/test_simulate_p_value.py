""" unit test the WeightedChoice class"""

from __future__ import division

import os
import unittest
import math

from denovonear.weights import WeightedChoice
from denovonear.transcript import Transcript
from denovonear.site_specific_rates import SiteRates
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.simulate import get_p_value

class TestGetPValuePy(unittest.TestCase):
    """ unit test the simulation of p-values
    """
    
    def setUp(self):
        
        self.transcript = self.construct_gene()
        self.rates = self.get_rates(self.transcript)
    
    def construct_gene(self):
        
        chrom = "1"
        name = "TEST"
        strand = "+"
        start = 0
        end = 70
        exons = [(5, 58)]
        cds = [(5, 58)]
        transcript = Transcript(name, chrom, start, end, strand)
        transcript.set_exons(exons, cds)
        transcript.set_cds(cds)
        
        cds = "ATGTGGGCTCCACCAGCAGCAATCATGGGATGGGCCCACCAAGAAGGTGGGTAA"
        gdna = "GGGGGATGTGGGCTCCACCAGCAGCAATCATGGGATGGGCCCACCAAGAAGGTGGGTAACCAGGCCCC"
        transcript.add_cds_sequence(cds)
        transcript.add_genomic_sequence(gdna)
        
        return transcript
    
    def get_rates(self, tx):
        # load the sequence contect mutation rates, then assess each site in the
        # CDS.
        mut_dict = load_mutation_rates()
        
        return SiteRates(tx, mut_dict)
    
    def test_get_p_value_single_de_novo(self):
        """ check that get_p_value() works correctly
        """
        
        # check that we don't assess transcripts with a single de novo
        iterations = 10000
        cq = 'missense'
        de_novos = [5]
        
        (obs, p_value) = get_p_value(self.transcript, self.rates, iterations, cq, de_novos)
        self.assertTrue(math.isnan(obs))
        self.assertTrue(math.isnan(p_value))
    
    def test_get_p_value(self):
        """
        """
        
        iterations = 10000
        cq = 'missense'
        de_novos = [5, 5]
        (obs, p_value) = get_p_value(self.transcript, self.rates, iterations, cq, de_novos)
        
        self.assertTrue(p_value < 0.04)
        self.assertEqual(obs, '0.0')
    
    def test_get_p_value_nonsignificant(self):
        """ check for de novos spread across the gene
        """
        
        iterations = 10000
        cq = 'missense'
        de_novos = [5, 58]
        
        (obs, p_value) = get_p_value(self.transcript, self.rates, iterations, cq, de_novos)
        
        self.assertTrue(p_value == 1.0)
        self.assertEqual(obs, '53.0')
    
    def test_get_p_value_lofs(self):
        """ check for loss of function rates
        """
    
        # construct two different tx objects, for the same gene, since otherwise
        # the p-values from different runs are the same
        tx1 = self.construct_gene()
        rates1 = self.get_rates(tx1)
        
        iterations = 10000
        cq = 'loss_of_function'
        de_novos = [5, 6]
        (obs_1, p_1) = get_p_value(tx1, rates1, iterations, cq, de_novos)
        
        tx2 = self.construct_gene()
        rates2 = self.get_rates(tx2)
        
        # make sure we can use the string 'lof' to get loss-of-function rates.
        # NOTE: due to randomly sampling, this will fail ~0.1% of the time,
        # purely by chance. If this fails, first try rerunning the tests.
        iterations = 10000
        cq = 'lof'
        (obs_2, p_2) = get_p_value(tx2, rates2, iterations, cq, de_novos)
        self.assertEqual(obs_1, obs_2)
        self.assertTrue(abs(p_1 - p_2) < 0.017)
