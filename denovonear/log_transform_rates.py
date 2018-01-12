
import math

def log_transform(rates):
    """ log transform a numeric value, unless it is zero, or negative
    """
    
    transformed = []
    for key in ['missense', 'nonsense', 'splice_lof', 'splice_region',
            'synonymous']:
        
        try:
            value = math.log10(rates[key])
        except ValueError:
            value = "NA"
        except KeyError:
            continue
        
        transformed.append(value)
    
    return '\t'.join(map(str, transformed))
