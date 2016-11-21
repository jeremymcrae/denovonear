""" load trincleotide mutation rates
"""

from pkg_resources import resource_filename

def load_mutation_rates(path=None):
    """ load sequence context-based mutation rates
    
    Args:
        path: path to table of sequence context-based mutation rates. If None,
            this defaults to per-trinucleotide rates provided by Kaitlin Samocha
            (Broad Institute).
    
    Returns:
        list of [initial, changed, rate] lists e.g. [['AGA', 'ATA', '5e-8']]
    """
    
    if path is None:
        path = resource_filename(__name__, "data/rates.txt")
    
    rates = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("from"): # ignore the header line
                continue
            
            line = [ x.encode('utf8') for x in line.strip().split() ]
            rates.append(line)
    
    return rates
