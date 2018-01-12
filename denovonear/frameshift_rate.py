
import math

def include_frameshift_rates(path):
    """ add per-gene frameshift mutation rates to the output file
    
    We make a crude estimate of the frameshift mutation rate by assessing the
    total frameshift burden (across all genes in a dataset) as 1.25 times the
    nonsense burden. Then apportion the frameshift burden proportional to each
    genes CDS length.
    
    This as as per Nature Genetics 46:944-950 (2014) doi:10.1038/ng.3050.
    Note that this is an approximation, and doesn't allow for divergence due to
    the base compostion of the CDS.
    
    Args:
        path: path to the output mutation rates (for the nonsense, missense etc)
    """
    
    with open(path) as handle:
        lines = [ x.strip().split('\t') for x in handle ]
    
    nonsense = sum([ 10**(float(x[4])) for x in lines[1:] if x[4] != 'NA' ])
    length = sum([ int(x[2]) for x in lines[1:] if x[2] != 'NA' ])
    
    # add the frameshift rates to each line in turn, while writing the output
    # back to the output path
    frameshift_sum = nonsense * 1.25
    with open(path, "w") as handle:
        for line in lines:
            if line[0] == "transcript_id":
                line.append("frameshift_rate")
            elif line[4] == 'NA':
                line.append('NA')
            else:
                # estimate the frameshift rate for the gene
                frameshift = (float(line[2])/length) * frameshift_sum
                frameshift = math.log10(frameshift)
                line.append(str(frameshift))
            
            line = "\t".join(line) +"\n"
            handle.write(line)
