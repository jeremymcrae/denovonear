""" function to translate DNA codons to single character amino acids
"""

from __future__ import division
from __future__ import print_function

import gzip
import os


def load_conservation_scores(conservation_folder, chrom, start, end):
    """ load conservation scores for a chromosome region
    """
    
    filename = "chr{0}.phyloP46way.wigFix.gz".format(chrom)
    path = os.path.join(conservation_folder, filename)
    
    scores = {}
    with gzip.open(path, "rt") as f:
        header = f.readline()
        header = header.strip().split()
        if header[0] != "fixedStep":
            raise ValueError("check wig format of file")
        
        wig_start = int(header[2].split("=")[1])
        wig_step = int(header[3].split("=")[1])
        
        pos = wig_start
        for line in f:
            pos += 1
            
            if start < pos < end:  
                scores[pos] = float(line.strip())
            elif pos > end: # break out once we're past the gene
                break
        
    return scores


