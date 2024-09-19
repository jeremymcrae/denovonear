""" class to analyse clustering of known de novos in genes according to their
distances apart within the gene, and compare that to simulated de novo events
within the same gene.
"""

from typing import Dict, List, Tuple

from gencodegenes.transcript import Transcript

from denovonear.structures import valid_structure
from denovonear.site_specific_rates import SiteRates
from denovonear.weights import (geomean,
                                geomean_double,
                                get_distances,
                                get_structure_distances,
                                simulate_clustering,
                                simulate_structure_clustering,
                                )

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
    sim_prob = simulate_clustering(weights, iterations, len(de_novos), observed)
    
    observed = "{0:0.1f}".format(observed)
    
    return (observed, sim_prob)

def get_structure_p_value(transcript: Transcript,
                          rates: SiteRates,
                          structure_coords: Dict[Tuple[str, int], Dict[str, float]],
                          iterations: int,
                          consequence: str,
                          de_novos: List[int]) -> Tuple[str, float]:
    """ find the probability of getting de novos with a mean conservation
    
    The probability is the number of simulations where the mean conservation
    between simulated de novos is less than the observed conservation.
    
    Args:
        transcript: Transcript object for the current gene.
        rates: SiteRates object, which contains WeightedChoice entries for
            different consequence categories.
        iterations: number of simulations to perform
        structure_coords: 3D coordinates for carbon atoms in amino acids along protein 
            sequence, indexed by (chain, residue number)
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
    
    if not valid_structure(structure_coords, transcript):
        return float('nan'), float('nan')
    
    coords = {k[1]: v for k, v in sorted(structure_coords.items())}
    cds_coords = [transcript.get_coding_distance(x)['pos'] // 3 for x in de_novos]
    cds_coords = [coords[x] for x in cds_coords]
    
    distances = get_structure_distances(cds_coords)
    observed = geomean_double(distances)
    
    coords = list(coords.values())
    # call a cython wrapped C++ library to handle the simulations
    sim_prob = simulate_structure_clustering(weights, coords, iterations, len(de_novos), observed)
    
    observed = "{0:0.1f}".format(observed)
    
    return (observed, sim_prob)
