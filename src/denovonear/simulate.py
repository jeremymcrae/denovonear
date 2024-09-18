""" class to analyse clustering of known de novos in genes according to their
distances apart within the gene, and compare that to simulated de novo events
within the same gene.
"""

import logging

from denovonear.weights import (geomean,
                                get_distances,
                                get_structure_distances,
                                analyse_de_novos,
                                analyse_structure_de_novos)

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
    
    cds_positions = [ transcript.get_coding_distance(x)['pos'] for x in de_novos ]
    distances = get_distances(cds_positions)
    observed = geomean(distances)
    
    # call a cython wrapped C++ library to handle the simulations
    sim_prob = analyse_de_novos(weights, iterations, len(de_novos), observed)
    
    observed = "{0:0.1f}".format(observed)
    
    return (observed, sim_prob)

def get_structure_p_value(transcript, rates, structure_coords, iterations, consequence, de_novos):
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
    
    # check there is only a single protein chain in the structure, fail if not
    n_chains = len(set(x[0] for x in structure_coords))
    if n_chains > 1:
        logging.warning(f'cannot get distances from structure with multiple chains')
        return float('nan'), float('nan')
    
    # make sure coords are linearly sorted without gaps, and that the protein
    # coords extend to the end of the transcript range
    residues = [x[1] for x in sorted(structure_coords)]
    if residues != list(range(1, len(structure_coords) + 1)):
        logging.warning(f'cannot get distances from structure missing residues')
        return float('nan'), float('nan')
    
    last_residue = transcript.get_coding_distance(transcript.get_cds_end())['pos'] // 3
    if abs(last_residue - residues[-1]) > 2:
        logging.warning(f"transcript length doesn't match structure length")
        return float('nan'), float('nan')
    
    coords = [v for k, v in sorted(structure_coords.items())]
    cds_coords = [ coords[transcript.get_coding_distance(x)['pos'] // 3] for x in de_novos ]
    
    distances = get_structure_distances(cds_coords)
    observed = geomean(distances)
    
    # call a cython wrapped C++ library to handle the simulations
    sim_prob = analyse_structure_de_novos(weights, coords, iterations, len(de_novos), observed)
    
    observed = "{0:0.1f}".format(observed)
    
    return (observed, sim_prob)
