""" unit test the WeightedChoice class"""

from __future__ import division

import unittest

from denovonear.weighted_choice import WeightedChoice

class TestWeightedChoicePy(unittest.TestCase):
    """ unit test the WeightedChoice class
    """
    
    def sample(self, chooser, max_iter):
        """ repeatedly sample a weighted choice object, up to max_iter
        """
        
        samples = []
        
        iteration = 0
        while iteration < max_iter:
            iteration += 1
            samples.append(chooser.choice())
        
        return samples
        
    def test___init__(self):
        """ check that __init__() sets the cumulative probability correctly
        """
        
        # check for two values
        choices = WeightedChoice([("a", 1), ("b", 5)])
        self.assertEqual(choices.cum_probs, [1, 6])
         
         # check for three values
        choices = WeightedChoice([("a", 1), ("b", 5), ("c", 10)])
        self.assertEqual(choices.cum_probs, [1, 6, 16])
        
        # check that it works for unsorted probabilities
        choices = WeightedChoice([("a", 1), ("b", 10), ("c", 5)])
        self.assertEqual(choices.cum_probs, [1, 11, 16])
        
        # check for very low values, with very high precision (but not
        # necessarily exactly equal)
        choices = WeightedChoice([("a", 5e-9), ("b", 1e-8), ("c", 1.000000000000005e-10)])
        self.assertAlmostEqual(choices.cum_probs[0], 5e-9, places=23)
        self.assertAlmostEqual(choices.cum_probs[1], 1.5e-8, places=23)
        self.assertAlmostEqual(choices.cum_probs[2], 1.51000000000000005e-8, places=23)
    
    def test_choice(self):
        """ test that choice() works correctly.
        
        Since WeightedChoice is a weighted random sampler, we can't rely on
        getting exact values out, so repeated samples are expected to obtain
        proportions of values equivalent to their weight value. The difference
        to the expected proportion minimises with larger sample sets, but at
        the cost of making the test hang for > 1 second for 1 million seconds,
        or > 10 s for 10 million samples.
        """
        
        choices = WeightedChoice([("a", 1), ("b", 5)])
        s = self.sample(choices, 2000000)
        self.assertAlmostEqual(s.count("a")/len(s), 0.1667, places=3)
        
        choices = WeightedChoice([("a", 1), ("b", 5), ("c", 4)])
        s = self.sample(choices, 100000)
        self.assertAlmostEqual(s.count("a")/len(s), 0.100, places=2)
        self.assertAlmostEqual(s.count("b")/len(s), 0.500, places=2)
        self.assertAlmostEqual(s.count("c")/len(s), 0.400, places=2)
        
        # very small numbers at the end still have expected proportions
        choices = WeightedChoice([("a", 1), ("b", 5), ("c", 0.0001)])
        s = self.sample(choices, 100000)
        self.assertAlmostEqual(s.count("c")/len(s), 0.0001, places=3)
        
        # very small numbers at the start still have expected proportions
        choices = WeightedChoice([("a", 0.0001), ("b", 1), ("c", 5)])
        s = self.sample(choices, 100000)
        self.assertAlmostEqual(s.count("a")/len(s), 0.0001, places=3)
        
        # check that the sampling works correctly at low weight values
        small = [x * 0.000000000001 for x in range(1000, 3000)]
        choices = WeightedChoice(list(zip(small, small)))
        s = self.sample(choices, 100000)
        self.assertAlmostEqual(s.count(small[0])/len(s), 0.0001, places=3)
        
        
    

# if __name__ == '__main__':
#     unittest.main()
