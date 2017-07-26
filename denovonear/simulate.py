""" class to analyse clustering of known de novos in genes according to their
distances apart within the gene, and compare that to simulated de novo events
within the same gene.
"""

from denovonear.weights import geomean, get_distances, analyse_de_novos

def get_p_value(transcript, rates, iterations, consequence, de_novos):
    """ find the probability of getting de novos with a mean conservation
    
    The probability is the number of simulations where the mean conservation
    between simulated de novos is less than the observed conservation.
    
    Args:
        transcript: Transcript object for the current gene.
        rates: SiteRates object, which contains WeightedChoice entries for
            different consequence categories.
        iterations: number of simulations to perform
        consequence: string to indicate the consequence type e.g. "missense, or
            "lof", "synonymous" etc. The full list is "missense", "nonsense",
            "synonymous", "lof", "loss_of_function", "splice_lof",
            "splice_region".
        de_novos: list of de novos within a gene
    
    Returns:
        tuple of mean proximity for the observed de novos and probability of
        obtaining a value less than or equal to the observed proximity from the
        null distribution.
    """
    
    if len(de_novos) < 2:
        return (float('nan'), float('nan'))
    
    rename = {"lof": "loss_of_function"}
    if consequence in rename:
        consequence = rename[consequence]
    
    weights = rates[consequence]
    
    cds_positions = [ transcript.chrom_pos_to_cds(x)['pos'] for x in de_novos ]
    distances = get_distances(cds_positions)
    observed = geomean(distances)
    
    # call a cython wrapped C++ library to handle the simulations
    sim_prob = analyse_de_novos(weights, iterations, len(de_novos), observed)
    
    observed = "{0:0.1f}".format(observed)
    
    return (observed, sim_prob)
