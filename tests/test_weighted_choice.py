""" unit test the WeightedChoice class"""

from __future__ import division

import unittest

from denovonear.weights import WeightedChoice

class TestWeightedChoicePy(unittest.TestCase):
    """ unit test the WeightedChoice class
    """
        
    def test___init__(self):
        """ check that __init__() initiates the object correctly
        """
        
        choices = WeightedChoice()
        
        # check that an object without any possible choices has a cumulative
        # sum of 0, but returns a choice of -1
        self.assertEqual(choices.get_summed_rate(), 0)
        self.assertEqual(choices.choice(), -1)
        
        # check that the type is set correctly
        self.assertEqual(type(choices), WeightedChoice)
        
    def test_add_choice(self):
        """ test that add_choice() works correctly
        """
        
        # check the cumulative sum while adding in new values
        choices = WeightedChoice()
        choices.add_choice(1, 1)
        self.assertEqual(choices.get_summed_rate(), 1)
        choices.add_choice(2, 5)
        self.assertEqual(choices.get_summed_rate(), 6)
        choices.add_choice(3, 10)
        self.assertEqual(choices.get_summed_rate(), 16)
        
        # check that it works for unsorted probabilities
        choices = WeightedChoice()
        choices.add_choice(1, 1)
        choices.add_choice(2, 10)
        choices.add_choice(3, 5)
        self.assertEqual(choices.get_summed_rate(), 16)
        
        # check for very low values, with very high precision (but not
        # necessarily exactly equal)
        choices = WeightedChoice()
        choices.add_choice(1, 5e-9)
        choices.add_choice(2, 1e-8)
        choices.add_choice(3, 1.000000000000005e-10)
        self.assertAlmostEqual(choices.get_summed_rate(), 1.51000000000000005e-8, places=23)
    
    def test_append(self):
        """ test that append() works correctly
        """
        
        # construct two objects
        a = WeightedChoice()
        a.add_choice(1, 0.5)
        
        b = WeightedChoice()
        b.add_choice(2, 1)
        
        # add one object to the other
        a.append(b)
        
        # check that the first object has changed correctly, but the other
        # remains unchanged
        self.assertEqual(a.get_summed_rate(), 1.5)
        self.assertEqual(b.get_summed_rate(), 1.0)


    
    def test_choice(self):
        """ test that choice() works correctly.
        
        Since WeightedChoice is a weighted random sampler, we can't rely on
        getting exact values out, so repeated samples are expected to obtain
        proportions of values equivalent to their weight value. The difference
        to the expected proportion minimises with larger sample sets, but at
        the cost of making the test hang for > 1 second for 1 million samples,
        or > 10 s for 10 million samples.
        """
        
        iterations = 1000000
        
        choices = WeightedChoice()
        choices.add_choice(1, 1)
        choices.add_choice(2, 5)
        s = [ choices.choice() for x in range(iterations) ]
        self.assertAlmostEqual(s.count(1)/len(s), 0.1667, places=2)
        
        # add another choice, then check that all of the choices have been
        # sampled at the expecetd proportions
        choices.add_choice(3, 4)
        s = [ choices.choice() for x in range(iterations) ]
        self.assertAlmostEqual(s.count(1)/len(s), 0.100, places=2)
        self.assertAlmostEqual(s.count(2)/len(s), 0.500, places=2)
        self.assertAlmostEqual(s.count(3)/len(s), 0.400, places=2)
        
        # check that all the choices have been made from the inserted values
        self.assertEqual(set(s), set([1, 2, 3]))
    
    def test_choice_small_numbers(self):
        """ test that choice() works correctly.
        """
        
        iterations = 1000000
        
        # very small numbers at the end still have expected proportions
        choices = WeightedChoice()
        choices.add_choice(1, 1)
        choices.add_choice(2, 5)
        choices.add_choice(3, 0.0001)
        s = [ choices.choice() for x in range(iterations) ]
        self.assertAlmostEqual(s.count(3)/len(s), 0.0001, places=3)
        
        # very small numbers at the start still have expected proportions
        choices = WeightedChoice()
        choices.add_choice(1, 0.0001)
        choices.add_choice(2, 1)
        choices.add_choice(3, 5)
        s = [ choices.choice() for x in range(iterations) ]
        self.assertAlmostEqual(s.count(1)/len(s), 0.0001, places=3)
        
        # check that the sampling works correctly at low weight values
        choices = WeightedChoice()
        
        numbers = range(1000, 3000)
        small = [ x * 0.000000000001 for x in numbers ]
        for (name, prob) in zip(numbers, small):
            choices.add_choice(name, prob)
        
        s = [ choices.choice() for x in range(iterations) ]
        self.assertAlmostEqual(s.count(numbers[0])/len(s), 0.0001, places=3)
    
    def test_choice_with_alleles(self):
        """ test that choice_with_alleles() works correctly.
        """
        
        # if you add a choice with alleles, then check that we get back alleles,
        # and that they are the same
        choices = WeightedChoice()
        choices.add_choice(1, 1, "A", "T")
        self.assertEqual(choices.choice_with_alleles(),
            {'alt': 'T', 'ref': 'A', 'pos': 1, 'offset': 0})
        self.assertEqual(choices.choice(), 1)
        
        # if you add choices without alleles, then default the alleles to "N"
        choices = WeightedChoice()
        choices.add_choice(1, 1)
        self.assertEqual(choices.choice_with_alleles(),
            {'alt': 'N', 'ref': 'N', 'pos': 1, 'offset': 0})
        
        # make sure you can't add multi-base alleles to the choices
        with self.assertRaises(TypeError):
            choices.add_choice(1, 1, "AA", "A")
            choices.add_choice(1, 1, "A", "AG")
        
        # make sure non-zero offsets are returned corectly
        choices = WeightedChoice()
        choices.add_choice(1, 1, "A", "T", 3)
        self.assertEqual(choices.choice_with_alleles(),
            {'alt': 'T', 'ref': 'A', 'pos': 1, 'offset': 3})
        self.assertEqual(choices.choice(), 1)
