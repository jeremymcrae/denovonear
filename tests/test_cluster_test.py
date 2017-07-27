""" unit test the  class"""

import os
import unittest
import tempfile
import shutil
import math

from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.cluster_test import fishers_method, cluster_de_novos, combine_p_values

class TestClusterTestPy(unittest.TestCase):
    """ unit test the genic test
    """
    
    @classmethod
    def setUpClass(cls):
        cls.temp_dir = tempfile.mkdtemp()
        cls.ensembl = EnsemblRequest(cls.temp_dir, 'grch37')
        cls.mut_dict = load_mutation_rates()
    
    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_dir)
    
    def test_cluster_de_novos(self):
        """ check that cluster_de_novos() runs correctly
        """
        
        symbol = 'PPP2R5D'
        de_novos = {'missense': [42975003, 42975003, 42975003, 42975013],
            'nonsense': []}
        
        p_values = cluster_de_novos(symbol, de_novos, 1000000,
            self.ensembl, self.mut_dict)
        
        self.assertAlmostEqual(p_values['miss_prob'], 3e-06, delta=3e-6)
        self.assertEqual(p_values['miss_dist'], '2.3')
        self.assertTrue(math.isnan(p_values['nons_prob']))
        self.assertEqual(p_values['nons_dist'], 'nan')
        
        # TODO: add a test with de novos that occur on different transcripts
    
    def test_fishers_method(self):
        """ check that fishers combined method calculates correctly
        """
        
        p_values = [0.01, 0.001]
        self.assertAlmostEqual(fishers_method(p_values), 0.0001251292546)
        
        # correct for single value. This is off by a minuscule fraction due to
        # float division
        p_values = [0.01]
        self.assertAlmostEqual(fishers_method(p_values), 0.01)
        
        # correct for list with NA value
        p_values = [0.01, float('nan')]
        self.assertAlmostEqual(fishers_method(p_values), 0.01)
        
        # correct for list with NA value
        p_values = [0.01, 0.001, float('nan')]
        self.assertAlmostEqual(fishers_method(p_values), 0.0001251292546)
        
        # correct without any p-values
        p_values = [float('nan'), float('nan')]
        self.assertTrue(math.isnan(fishers_method(p_values)))
        
        # raise error for list with zero
        with self.assertRaises(ValueError):
            p_values = [0, 0.01]
            fishers_method(p_values)
    
    def test_combine_p_values(self):
        """ check that combine_p_values() works correctly
        """
        
        probs = {"miss_prob": [0.01, 0.001], "nons_prob": [0.01, 0.01]}
        data = combine_p_values(probs)
        
        self.assertAlmostEqual(data['miss_prob'], 0.0001251292546)
        self.assertAlmostEqual(data['nons_prob'], 0.0010210340371)
