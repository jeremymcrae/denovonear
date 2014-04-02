""" class to analyse clustering of known de novos in genes according to their 
distances apart within the gene, and compare that to simulated de novo events 
within the same gene.
"""

from __future__ import division
from __future__ import print_function

import bisect
import itertools
import math
import operator
import random


from analyse_de_novos import AnalyseDeNovos

class AnalyseDeNovoConservation(AnalyseDeNovos):
    """ class to analyse conservation of de novo events via site specific 
    mutation rates
    """
    
    def get_score(self, positions):
        """ gets the mean conservation of two or more CDS positions
        
        Args:
            positions: list of numbers
        
        Returns:
            provides the mean conservation at the gene positions
        """
        
        assert len(positions) > 1
        
        scores = []
        for position in positions:
            score = self.transcript.get_conservation_score(position)
            scores.append(score)
        
        # get the mean score,  but make the score negative, since the simulation
        # looks for a minimal value
        return -(sum(scores)/len(scores))
    
    def analyse_de_novos(self, de_novos, weights, iterations):
        """
        """
        
        if len(de_novos) < 2:
            return ("NA", "NA")
        
        minimum_prob = 1/(1 + iterations)
        sim_prob = minimum_prob
        dist = []
        
        cds_positions = self.convert_de_novos_to_cds_positions(de_novos)
        observed_distances = self.get_log_distances(cds_positions)
        mean_distance = 10 ** (sum(observed_distances)/len(observed_distances))
        
        observed_value = self.get_score(cds_positions)
        
        # if the p-value that we obtain is right at the minimum edge of the 
        # simulated distribution, increase the number of iterations until the
        # p-value is no longer at the very edge (or we reach 100 million 
        # iterations).
        while iterations < 100000000 and sim_prob == minimum_prob:
            minimum_prob = 1/(1 + iterations)
            
            dist = self.simulate_distribution(de_novos, mean_distance, weights, dist, len(de_novos), iterations)
            pos = bisect.bisect_right(dist, observed_value)
            sim_prob = (1 + pos)/(1 + len(dist))
            
            iterations += 1000000 # for if we need to run more iterations
        
        # output = open("/nfs/users/nfs_j/jm33/apps/mutation_rates/data/distribution.txt", "w")
        # for val in dist:
        #     output.write(str(val) + "\n")
        
        if type(observed_value) != "str":
            observed_value = "{0:0.1f}".format(observed_value)
        
        return (observed_value, sim_prob)
        
    def simulate_distribution(self, de_novos, observed, weights, dist=[], sample_n=2, max_iter=100):
        """ creates a distribution of mutation scores in a single gene
        
        Args:
            weights: WeightedChoice object
            dist: current list of simulated scores
            sample_n: number of de novo mutations to sample
            max_iter: number of iterations/simulations to run
        """
        
        iteration = len(dist)
        while iteration < max_iter:
            iteration += 1
            
            positions = self.get_positions(de_novos, observed, weights)
            
            # the following line is class specific - can do distance clustering,
            # conservation scores
            value = self.get_score(positions)
            dist.append(value)
        
        dist.sort()
        
        return dist
    
    def get_positions(self, de_novos, observed_value, weights):
        """ get positions that match the observed clustering distance
        """
        
        weighted_sites = weights.get_positions()
        
        cds_start = 1
        cds_end = len(self.transcript.cds_sequence)
        
        positions = []
        sample = random.choice(list(weighted_sites))
        positions.append(sample)
        
        while len(positions) < len(de_novos):
            allowed = self.find_positions(positions, len(de_novos), observed_value, cds_start, cds_end)
            allowed = allowed & weighted_sites
            
            # sometimes we cannot find suitable sites (because we ignore ranges 
            # where both end points are unsuitable, and sometimes the end points 
            # are unsuitable, but some middle points are suitable).
            # The quickest fix is to resample from the beginning
            while len(allowed) == 0:
                return self.get_positions(de_novos, observed_value, weights)
            
            if len(allowed) > 1:
                sample = weights.choose_amongst_positions(allowed)
            else:  
                sample = list(allowed)[0]
            positions.append(sample)
        
        return positions

    def get_log_distances(self, sites):
        """ gets a list of log10 distances between sites
        """
        
        pos_pairs = itertools.combinations(sites, 2)
        distances = []
        for pos_1, pos_2 in pos_pairs:
            distance = abs(pos_1 - pos_2)
            if distance == 0:
                distance = 1
            distances.append(math.log10(distance))
        
        return distances
    
    def get_segments(self, sites, start, end):
        """ split a list of sites into the segments that separate those sites
        """
        
        sites = sorted(set(sites))
        segments = []
        for pos in range(len(sites) + 1):
            if pos == 0:
                segments.append((start, sites[0], "left"))
            elif pos == len(sites):
                segments.append((sites[-1] + 1, end, "right"))
            else:
                midpoint = (sites[pos - 1] + sites[pos])//2
                segments.append((sites[pos - 1] + 1, midpoint, "right"))
                segments.append((midpoint + 1, sites[pos], "left"))
        
        return segments
    
    def find_midpoint(self, high, low, side, sites, dist, log_dist):
        """ bisects a range to find a sequence of positions
        """
        
        midpoint = (high + low)/2
        while abs(high - low) > 0.4:
            midpoint = (high + low)/2
            
            passes = self.check_position(midpoint, sites, dist, log_dist, False)
            if passes:
                if side == "left":
                    high = midpoint
                else:
                    low = midpoint
            else:
                if side == "left":
                    low = midpoint
                else:
                    high = midpoint
        
        midpoint = int(round(midpoint))
        
        return midpoint
    
    def find_positions(self, sites, de_novos, observed, start, end):
        """ find suitable positions to match the clustering distance
        """
        
        # avoid taking the log of 0
        if observed == 0:
            observed = 1
        
        log_dist = math.log10(observed)
        repeat_site = random.randint(0, len(sites) - 1)
        
        current_n = len(sites)
        extended_sites = sites + [sites[repeat_site]] * ((de_novos - len(sites)) - 1)
        dist = self.get_log_distances(extended_sites)
        
        segments = self.get_segments(extended_sites, start, end)
        
        ranges = set([])
        for segment in segments:
            low = segment[0]
            high = segment[1]
            side = segment[2]
            
            if not self.check_position(low, extended_sites, dist, log_dist, False) and \
                   not self.check_position(high, extended_sites, dist, log_dist, False):
                continue
            
            midpoint = self.find_midpoint(high, low, side, extended_sites, dist, log_dist)
            
            if high == low:
                subrange = range(low, high + 1)
            if side == "left":
                subrange = range(midpoint, segment[1] + 1)
            else:
                subrange = range(segment[0], midpoint + 1)
            
            ranges = ranges | set(subrange)
        
        # if we are selecting the final position, we want to find the exact site 
        # that will take us to the correct average distance, so look through the 
        # suitable range positions for sites that give the exact average (or at 
        # least as close as we will tolerate)
        if de_novos - current_n == 1:
            single_positions = set([])
            for pos in ranges:
                if self.check_position(pos, extended_sites, dist, log_dist, exact=True):
                    single_positions.add(pos)
            ranges = single_positions
        
        return ranges
    
    def get_exact_sites(self):
        pass
    
    def check_position(self, new_pos, sites, dist, log_dist, exact=False):
        """ checks if a position can enable the desired clustering distance
        """
        
        new_dist = dist[:]
        
        # get the remaining distance pairs
        for pos in sites:
            difference = abs(new_pos - pos)
            if difference == 0:
                difference = 1
            new_dist.append(math.log10(difference))
        
        # mostly we only want to check if the position won't exceed the desired
        # distance
        if not exact:
            if sum(new_dist)/len(new_dist) <= log_dist:
                return True
        else:
            # allow us to only use sites if they get us to the exact desired
            # value (can tolerate up to 3 bp difference)
            required = 10 ** log_dist
            current = 10 ** (sum(new_dist)/len(new_dist))
            if abs(current - required) <= 2:
                return True
        
        return False
    
