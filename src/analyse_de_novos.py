""" class to analyse clustering of known de novos in genes according to their
distances apart within the gene, and compare that to simulated de novo events
within the same gene.
"""

from __future__ import division
from __future__ import print_function

import bisect
import math
import ctypes
import os
import sys
import heapq

IS_PYTHON2 = sys.version_info[0] == 2
IS_PYTHON3 = sys.version_info[0] == 3

BUILD_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "build")

class AnalyseDeNovos(object):
    """ class to analyse clustering of de novo events via site specific
    mutation rates
    """
    
    def __init__(self, transcript, site_weights, iterations):
        """ initialise the class
        """
        
        # define the c library to use
        if IS_PYTHON2:
            lib_path = os.path.join(BUILD_DIR, "lib.linux-x86_64-2.7", "libsimulatedenovo.so")
        elif IS_PYTHON3:
            lib_path = os.path.join(BUILD_DIR, "lib.linux-x86_64-3.3", "libsimulatedenovo.cpython-33m.so")
        
        if os.path.exists(lib_path):
            self.lib = ctypes.cdll.LoadLibrary(lib_path)
            # make sure we set the return type
            self.lib.c_analyse_de_novos.restype = ctypes.c_double
            self.analyse_de_novos = self.c_analyse_de_novos
        
        self.transcript = transcript
        self.site_weights = site_weights
        self.max_iter = iterations
    
    def analyse_missense(self, de_novo_events):
        """ analyse clustering of missense de novos
        """
        
        weights = self.site_weights.get_missense_rates_for_gene()
        return self.analyse_de_novos(de_novo_events, weights, self.max_iter)
    
    def analyse_nonsense(self, de_novo_events):
        """ analyse clustering of nonsense de novos
        """
        
        weights = self.site_weights.get_nonsense_rates_for_gene()
        return self.analyse_de_novos(de_novo_events, weights, self.max_iter)
    
    def analyse_functional(self, de_novo_events):
        """ analyse clustering of functional (missense and nonsense) de novos
        """
        
        weights = self.site_weights.get_functional_rates_for_gene()
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
        self.dist = []
        
        cds_positions = self.convert_de_novos_to_cds_positions(de_novos)
        observed_value = self.get_score(cds_positions)
        
        # if the p-value that we obtain is right at the minimum edge of the 
        # simulated distribution, increase the number of iterations until the
        # p-value is no longer at the very edge (or we reach 100 million 
        # iterations).
        while iterations < 100000000 and sim_prob == minimum_prob:
            minimum_prob = 1/(1 + iterations)
            
            self.simulate_distribution(weights, len(de_novos), iterations)
            pos = bisect.bisect_right(self.dist, observed_value)
            sim_prob = (1 + pos)/(1 + len(self.dist))
            
            # assess whether the P-value could never fall below 0.1, and cut 
            # out after a smaller number of iterations, in order to minimise 
            # run time. Figure out the lower bound of the confidence interval
            # for the current simulated P value.
            z = 10
            alpha = 0.1
            delta = (z * math.sqrt((sim_prob * (1 - sim_prob))))/iterations
            lower_bound = sim_prob - delta
            
            # if the lower bound of the confidence interval exceeds 0.1, then we
            # can be sure it's not going to ever get lower than 0.05.
            if lower_bound > alpha:
                break
            
            iterations += 100000 # for if we need to run more iterations
        
        # output = open("/nfs/users/nfs_j/jm33/apps/mutation_rates/data/distribution.txt", "w")
        # for val in dist:
        #     output.write(str(val) + "\n")
        
        if type(observed_value) != "str":
            observed_value = "{0:0.1f}".format(observed_value)
        
        return (observed_value, sim_prob)
    
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
        observed_value = self.get_score(cds_positions)
        
        sites = []
        probs = []
        for (pos, rate) in weights.choices:
            sites.append(pos)
            probs.append(rate)
        
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
        sim_prob = self.lib.c_analyse_de_novos(c_sites, c_probs, length, iterations, de_novo_count, c_observed_value)
        
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
        
        # # occasionally we want to find the dispersion of sampled sites in a gene
        # # so we write the sampled sites to a file for analysis
        # src_dir = os.path.dirname(__file__)
        # cluster_dir = os.path.dirname(src_dir)
        # path = os.path.join(cluster_dir, "data", "sampled_sites.weighted.txt")
        # output = open(path, "w")
        # chrom_positions = []
        
        iteration = len(self.dist)
        temp_dist = []
        while iteration < max_iter:
            iteration += 1
            
            positions = []
            while len(positions) < sample_n:
                site = weights.choice()
                positions.append(site)
                # chr_pos = self.transcript.get_position_on_chrom(site)
                # chrom_positions.append(str(chr_pos))
            
            # the following line is class specific - can do distance clustering,
            # conservation scores
            value = self.get_score(positions)
            temp_dist.append(value)
        
        # output.write("\n".join(chrom_positions))
        # output.close()
        
        # sort the current distances, then use a
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
    
    def geomean(self, values):
        """ get the geometric mean of a list of values
        """
        
        # get the geometric mean, but be careful around values of 0, since
        # without correction, the mean distance would be zero
        if 0 in values:
            # allow for 0s in a geometric mean by shifting everything up one, 
            # then dropping the mean by one at the end
            values = [math.log10(x + 1) for x in values]
            logmean = sum(values)/len(values)
            mean = (10 ** logmean) - 1
        else:
            values = [math.log10(x) for x in values]
            logmean = sum(values)/len(values)
            mean = 10 ** logmean
        
        return mean


