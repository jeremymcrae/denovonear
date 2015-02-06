""" load trincleotide mutation rates
"""

from __future__ import print_function

def load_trincleotide_mutation_rates(filename, indel_only=False):
    """ load mutation rates provided by Kaitlin Samocha (Broad Institute).
    """
    
    mut_dict = {}
    with open(filename) as f:
        for line in f:
            if line.startswith("from"): # ignore the header line
                continue
            
            line = line.strip().split()
            initial = line[0]
            changed = line[1]
            rate = float(line[2])
            
            if indel_only:
                rate = 1e-8
            
            if initial not in mut_dict:
                mut_dict[initial] = {}
            
            mut_dict[initial][changed] = rate
    
    return mut_dict
