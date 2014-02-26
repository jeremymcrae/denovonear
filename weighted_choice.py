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
        
        if pos == len(self.cum_probs) - 1:
            return self.choices[pos][0]
        
        # since we have chosen the position immediately prior to the random
        # probability, figure out whether the value before or after is closer
        left_diff = abs(r - self.cum_probs[pos])
        right_diff = abs(r - self.cum_probs[pos + 1])
        
        if left_diff > right_diff:
            return self.choices[pos][0]
        else:
            return self.choices[pos + 1][0]
