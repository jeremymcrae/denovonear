""" class to analyse clustering of known de novos in genes according to their
distances apart within the gene, and compare that to simulated de novo events
within the same gene.
"""

from __future__ import division, print_function

import os
import sys
import ctypes
import glob

from denovonear.geometric_mean import geomean

class ClusterDeNovos(object):
    """ class to analyse clustering of de novos via site specific mutation rates
    """
    
    def __init__(self, transcript, rates, iterations):
        """ initialise the class
        
        Args:
            transcript: Transcript Object for the current gene.
            rates: SiteRates object, which contains WeightedChoice entries for
                different consequence categories.
            iterations: number of simulations to perform
        """
        
        # define the c library to use
        for path in sys.path:
            if os.path.isdir(path):
                files = glob.glob(os.path.join(path, "*simulatedenovo*so"))
                if len(files) > 0:
                    lib_path = files[0]
                    self.lib = ctypes.CDLL(lib_path)
                    
                    # make sure we set the return type
                    self.lib.analyse_de_novos.restype = ctypes.c_double
                    break
        
        self.transcript = transcript
        self.rates = rates
        self.iterations = iterations
        self.dist = []
    
    def analyse_de_novos(self, consequence, de_novos):
        """ find the probability of getting de novos with a mean conservation
        
        The probability is the number of simulations where the mean conservation
        between simulated de novos is less than the observed conservation.
        
        Args:
            consequence: string to indicate the consequence type e.g. "missense, or
                "lof", "synonymous" etc. The full list is "missense", "nonsense",
                "synonymous", "lof", "loss_of_function", "splice_lof",
                "splice_region".
            de_novos: list of de novos within a gene
        
        Returns:
            mean conservation for the observed de novos and probability of
            obtaining a mean conservation less than the observed conservation
        """
        
        if len(de_novos) < 2:
            return ("NA", "NA")
        
        rename = {"lof": "loss_of_function"}
        if consequence in rename:
            consequence = rename[consequence]
        
        weights = self.rates[consequence]
        
        cds_positions = self.convert_de_novos_to_cds_positions(de_novos)
        observed_value = geomean(cds_positions)
        
        # convert the number of iterations and de novos to ctypes
        iterations = ctypes.c_int(self.iterations)
        de_novo_count = ctypes.c_int(len(de_novos))
        c_observed_value = ctypes.c_double(observed_value)
        
        # call a C++ library to handle the simulations
        sim_prob = self.lib.c_analyse_de_novos(weights, iterations, de_novo_count,
            c_observed_value)
        
        if type(observed_value) != "str":
            observed_value = "{0:0.1f}".format(observed_value)
        
        return (observed_value, sim_prob)
    
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
