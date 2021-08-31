""" unit test the  class"""

import os
import unittest
import asyncio
import math

from denovonear.load_gene import load_gene
from denovonear.rate_limiter import RateLimiter
from denovonear.cluster_test import fishers_method, cluster_de_novos

async def call(func, *args, **kwargs):
    ''' call ensembl rest API function
    '''
    async with RateLimiter(15) as ensembl:
        return await func(*args, ensembl=ensembl, **kwargs)

def _run(func, *args, **kwargs):
    return asyncio.get_event_loop().run_until_complete(call(func, *args, **kwargs))

class TestClusterTestPy(unittest.TestCase):
    """ unit test the genic test
    """
    
    def test_cluster_de_novos(self):
        """ check that cluster_de_novos() runs correctly
        """
        
        symbol = 'PPP2R5D'
        de_novos = {'missense': [42975003, 42975003, 42975003, 42975013],
            'nonsense': []}
        
        gene = _run(load_gene, gene_id=symbol)
        p_values = cluster_de_novos(symbol, de_novos, gene, iterations=1000000)
        
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
