""" load trincleotide mutation rates
"""

def load_mutation_rates(filename):
    """ load mutation rates provided by Kaitlin Samocha (Broad Institute).
    """
    
    mut_dict = []
    with open(filename) as handle:
        for line in handle:
            if line.startswith("from"): # ignore the header line
                continue
            
            line = line.strip().split()
            initial = line[0]
            changed = line[1]
            rate = line[2]
            
            # if initial not in mut_dict:
            #     mut_dict[initial] = {}
            
            mut_dict.append([initial, changed, rate])
    
    return mut_dict
