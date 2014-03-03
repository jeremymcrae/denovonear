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
        pos = 0
        for line in f:
            pos += 1
            
            # occasionally the file jumps the position to a new site, not just
            # at the start of the file
            if line.startswith("fixed"):
                line = line.strip().split()
                pos = int(line[2].split("=")[1])
                continue
            
            if start < pos < end:  
                scores[pos] = float(line.strip())
            elif pos > end: # break out once we're past the gene
                break
        
    return scores


