""" class to analyse clustering of known de novos in genes according to their
distances apart within the gene, and compare that to simulated de novo events
within the same gene.
"""

from __future__ import division, print_function

import os
import sys
import bisect
import math
import ctypes
import glob
import heapq

from denovonear.geometric_mean import geomean

class ClusterDeNovos(object):
    """ class to analyse clustering of de novos via site specific mutation rates
    """
    
    def __init__(self, transcript, site_weights, iterations):
        """ initialise the class
        """
        
        # define the c library to use
        for path in sys.path:
            if os.path.isdir(path):
                files = glob.glob(os.path.join(path, "*simulatedenovo*so"))
                if len(files) > 0:
                    lib_path = files[0]
                    self.lib = ctypes.CDLL(lib_path)
                    
                    # make sure we set the return type
                    self.lib.c_analyse_de_novos.restype = ctypes.c_double
                    self.analyse_de_novos = self.c_analyse_de_novos
                    break
        
        self.transcript = transcript
        self.site_weights = site_weights
        self.max_iter = iterations
        self.dist = []
    
    def analyse_consequence(self, consequence):
        """ analyse clustering of a particulare functional consequence
        
        Args:
            category: string to indicate the consequence type e.g. "missense, or
                "lof", "synonymous" etc. The full list is "missense", "nonsense",
                 "functional", "synonymous", "lof", "loss_of_function",
                 "splice_lof", "splice_region", "missense_and_splice_region".
        
        Returns:
            p-value for getting the observing clustering within the functional
            category.
        """
        
        weights = self.site_weights.get_cds_rates(consequence)
        return self.analyse_de_novos(de_novo_events, weights, self.max_iter)
    
    def analyse_de_novos(self, de_novos, weights, iterations):
        """ find the probability of getting de novos with a mean conservation
        
        The probability is the number of simulations where the mean conservation
        between simulated de novos is less than the observed conservation.
        
        Args:
            de_novos: list of de novos within a gene
            weights: WeightedChoice object to randomly choose positions within
                a gene using site specific mutation rates.
            iterations: the (minimum) number of perumtations to run.
        
        Returns:
            mean conservation for the observed de novos and probability of
            obtaining a mean conservation less than the observed conservation
        """
        
        if len(de_novos) < 2:
            return ("NA", "NA")
        
        minimum_prob = 1/(1 + iterations)
        sim_prob = minimum_prob
        
        cds_positions = self.convert_de_novos_to_cds_positions(de_novos)
        observed_value = geomean(cds_positions)
        
        # if the p-value that we obtain is right at the minimum edge of the
        # simulated distribution, increase the number of iterations until the
        # p-value is no longer at the very edge (or we reach 100 million
        # iterations).
        while iterations < 100000000 and sim_prob == minimum_prob:
            minimum_prob = 1/(1 + iterations)
            
            self.simulate_distribution(weights, len(de_novos), iterations)
            pos = bisect.bisect_right(self.dist, observed_value)
            sim_prob = (1 + pos)/(1 + len(self.dist))
            
            # halt permutations if the P value could never be significant
            z = 10
            alpha = 0.1
            if self.adaptive_permutation(sim_prob, iterations, z, alpha):
                break
            
            iterations += 1000000 # for if we need to run more iterations
        
        if type(observed_value) != "str":
            observed_value = "{0:0.1f}".format(observed_value)
        
        return (observed_value, sim_prob)
    
    def adaptive_permutation(self, p_val, iterations, z_score=10, alpha=0.1):
        """ halt permutations if the P value could never be significant
        
        Args:
            p_val: current simulated P value
            iterations: iterations run in order to obtain the simulated P value
            z_score:
            alpha: threshold
        
        Returns:
            True/False for whether to halt the permuations
        """
        
        # assess whether the P-value could never fall below 0.1, and halt after
        # a smaller number of iterations, in order to minimise run time. Figure
        # out the lower bound of the confidence interval for the current
        # simulated P value.
        delta = (z_score * math.sqrt((p_val * (1 - p_val))))/iterations
        lower_bound = p_val - delta
        
        # if the lower bound of the confidence interval exceeds our alpha, then
        # we can be sure it's not going to ever get lower than 0.05.
        return lower_bound > alpha
    
    def c_analyse_de_novos(self, de_novos, weights, iterations):
        """ find the probability of getting de novos with a mean conservation
        
        The probability is the number of simulations where the mean conservation
        between simulated de novos is less than the observed conservation.
        
        Args:
            de_novos: list of de novos within a gene
            weights: WeightedChoice object to randomly choose positions within
                a gene using site specific mutation rates.
            iterations: the (minimum) number of perumtations to run.
        
        Returns:
            mean conservation for the observed de novos and probability of
            obtaining a mean conservation less than the observed conservation
        """
        
        if len(de_novos) < 2:
            return ("NA", "NA")
        
        cds_positions = self.convert_de_novos_to_cds_positions(de_novos)
        observed_value = geomean(cds_positions)
        
        sites, probs = zip(*weights.choices)
        
        # convert the sampler sites to c types
        c_sites = ctypes.c_int * len(sites)
        c_sites = c_sites(*sites)
        
        # convert the mutation rates associated with the sites into c types
        c_probs = ctypes.c_double * len(sites)
        c_probs = c_probs(*probs)
        
        # convert the number of iterations and de novos to ctypes
        length = ctypes.c_int(len(sites))
        iterations = ctypes.c_int(iterations)
        de_novo_count = ctypes.c_int(len(de_novos))
        c_observed_value = ctypes.c_double(observed_value)
        
        # call a C++ library to handle the simulations
        sim_prob = self.lib.c_analyse_de_novos(c_sites, c_probs, length,
            iterations, de_novo_count, c_observed_value)
        
        if type(observed_value) != "str":
            observed_value = "{0:0.1f}".format(observed_value)
        
        return (observed_value, sim_prob)
    
    def simulate_distribution(self, weights, sample_n=2, max_iter=100):
        """ creates a distribution of mutation scores in a single gene
        
        Args:
            weights: WeightedChoice object
            sample_n: number of de novo mutations to sample
            max_iter: number of iterations/simulations to run
        """
        
        iteration = len(self.dist)
        temp_dist = []
        while iteration < max_iter:
            iteration += 1
            
            positions = []
            while len(positions) < sample_n:
                site = weights.choice()
                positions.append(site)
            
            value = geomean(positions)
            temp_dist.append(value)
        
        # sort the current distances, then use a sort function that is fast
        # when merging two sorted lists
        temp_dist = sorted(temp_dist)
        self.dist = list(heapq.merge(self.dist, temp_dist))
    
    def convert_de_novos_to_cds_positions(self, de_novos):
        """ convert cds positions for de novo events into cds positions
        
        Args:
            de_novos: list of chrom bp positions within the transcript
        
        Returns:
            list of positions converted to CDS positions within the transcript
        """
        
        cds_positions = []
        for pos in de_novos:
            dist = self.transcript.convert_chr_pos_to_cds_positions(pos)
            cds_positions.append(dist)
        
        return cds_positions
