""" get a weighted random selection, from:
http://stackoverflow.com/questions/3679694/a-weighted-version-of-random-choice
"""

import random
import bisect

class WeightedChoice(object):
    """ class for weighted choices, rather than a function, so we don't have to 
    keep resumming the weights.
    """
    
    def __init__(self, choices):
        """ set up a list of cumulative probabilities
        """
        self.choices = choices
        
        # make a list of the cumulative probabilities
        cum_prob = 0
        self.cum_probs = []
        for pos in range(len(self.choices)):
            cum_prob += self.choices[pos][1]
            self.cum_probs.append(cum_prob)
        
    def choice(self):
        """ chooses a random position using a set of probability weights
        """
        
        # figure out where in the list a random probability would fall
        r = random.uniform(0, self.cum_probs[-1])
        pos = bisect.bisect_left(self.cum_probs, r)
        
        return self.choices[pos][0]
