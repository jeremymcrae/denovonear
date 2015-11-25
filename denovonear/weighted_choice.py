""" get a weighted random selection, from:
http://stackoverflow.com/questions/3679694/a-weighted-version-of-random-choice
"""

import random
import bisect

class WeightedChoice(object):
    """ class for weighted choices, rather than a function, so we don't have to
    keep resumming the weights.
    """
    
    def __init__(self, choices=None):
        """ set up a list of cumulative probabilities
        
        Args:
            choices: list of (name, probability) tuples, or None
        """
        self.choices = []
        self.cum_probs = []
        
        if choices is not None:
            for pair in choices:
                (name, prob) = pair
                self.add_choice(name, prob)
        
    def add_choice(self, name, probability):
        """ add another possible choice for selection
        
        Args:
            name: an ID, for example, CDS position, so we know each choice.
            probability: probability of selecting this choice.
        """
        
        self.choices.append((name, probability))
        
        cum_prob = self.get_summed_rate() + probability
        self.cum_probs.append(cum_prob)
    
    def get_summed_rate(self):
        """ return the cumulative probability for the class
        """
        
        if len(self.cum_probs) == 0:
            return 0
        else:
            return self.cum_probs[-1]
    
    def choice(self):
        """ chooses a random element using a set of probability weights
        
        Returns:
            the name of the randomly selected element (e.g. position)
        """
        
        # figure out where in the list a random probability would fall
        r = random.uniform(0, self.cum_probs[-1])
        pos = bisect.bisect_left(self.cum_probs, r)
        
        return self.choices[pos][0]
