
from math import log, isnan
from typing import Dict, List, Tuple, Union

from scipy.stats import chi2

from gencodegenes import Gene
from denovonear.load_gene import get_de_novos_in_transcript, minimise_transcripts
from denovonear.site_specific_rates import SiteRates
from denovonear.simulate import get_p_value, get_structure_p_value

def fishers_method(values):
    """ function to combine p values, using Fisher's method
    
    We occasionally have multiple P values for a mutation type, obtained from
    alternative transcripts for the gene. If we have only one P value for the
    gene for the mutation type, we just use that value, if we don't have any
    data, we use "NA", otherwise combine p-values using Fisher's method.
    
    Args:
        x: list of P-values for a gene
    
    Returns:
        combined P-value
    """
    values = [ x for x in values if not isnan(x) ]
    # use Fisher's combined method to estimate the P value from multiple
    # P-values. The chi square statistic is -2*sum(ln(P-values))
    return chi2.sf(-2 * sum(map(log, values)), 2 * len(values))

def cluster_de_novos(de_novos, gene, mut_dict, iterations=1000000):
    """ analysis proximity cluster of de novos in a single gene
    
    Args:
        de_novos: dictionary of de novo positions for the HGNC gene,
            indexed by functional type
        gene: gencodegenes.Gene object
        mut_dict: dictionary of mutation rates, indexed by trinuclotide sequence, 
        iterations: number of simulations to run
    
    Returns:
        a dictionary containing P values, and distances for missense, nonsense,
        and synonymous de novos events. Missing data is represented by "NA".
    """
    
    missense = de_novos["missense"]
    nonsense = de_novos["nonsense"]
    
    if len(gene.transcripts) == 0:
        nan = float('nan')
        return {'miss_dist': nan, 'miss_prob': nan, 'nons_prob': nan, 'nons_dist': nan}
    
    # load the set of transcripts that are the  minimum set of transcripts
    # required to contain all the de novos, unless we can't find any coding
    # transcripts that contain the de novos.
    transcripts = gene.transcripts
    minimized = minimise_transcripts(transcripts, missense + nonsense)
    transcripts = [x for x in transcripts if x.name in minimized]
    
    if len(transcripts) == 0:
        nan = float('nan')
        return {'miss_dist': nan, 'miss_prob': nan, 'nons_prob': nan, 'nons_dist': nan}
    
    probs = {"miss_prob": [], "nons_prob": []}
    dists = {"miss_dist": [], "nons_dist": []}
    
    for transcript in transcripts:
        
        missense_events = get_de_novos_in_transcript(transcript, missense)
        nonsense_events = get_de_novos_in_transcript(transcript, nonsense)
        
        # The rates either are passed in as lists of lists of sequence-context 
        # based rates which can be applied to a given sequence, or as 
        # dictionaries of rates per position/allele from genome based reference 
        # sources. If the latter, then we have to load the rates within the 
        # transcript region
        _rates = mut_dict
        if not isinstance(mut_dict, list):
            _rates = mut_dict(transcript.chrom, transcript.start, transcript.end)
        
        rates = SiteRates(transcript, _rates)
        
        (miss_dist, miss_prob) = get_p_value(transcript, rates, iterations, "missense", missense_events)
        (nons_dist, nons_prob) = get_p_value(transcript, rates, iterations, "lof", nonsense_events)
        
        dists["miss_dist"].append(miss_dist)
        dists["nons_dist"].append(nons_dist)
        probs["miss_prob"].append(miss_prob)
        probs["nons_prob"].append(nons_prob)
        
        # remove the de novos analysed in the current transcript, so that
        # analysis of subsequent transcripts uses independent events. NOTE THAT
        # THIS MIGHT MISS SOME CLUSTERING ACROSS MUTUALLY EXCLUSIVE TRANSCRIPTS
        # IF THE DE NOVO EVENTS ARE NEAR THE TRANSCRIPT DIVERGENCE.
        missense = [x for x in missense if x not in missense_events]
        nonsense = [x for x in nonsense if x not in  nonsense_events]
        
    for key in dists:
        dists[key] = ",".join([ str(x) for x in dists[key] ])
    
    probs = {k: fishers_method(probs[k]) for k in probs}
    probs.update(dists)
    
    return probs

def structure_cluster_de_novos(de_novos: Dict[str, List[int]],
                               structure: Dict[Tuple[str, int], Dict[str, float]],
                               gene: Gene,
                               mut_dict: Union[List[List[str]], Dict[int, Dict[str, float]]],
                               iterations=1000000,
                               ) -> Dict[str, float]:
    """ analysis proximity cluster of de novos in a single gene
    
    Args:
        de_novos: dictionary of de novo positions for the HGNC gene,
            indexed by functional type
        structure: 3D coordinates for carbon atoms in amino acids along protein 
            sequence, indexed by (chain, residue number)
        gene: gencodegenes.Gene object
        mut_dict: dictionary of mutation rates, indexed by trinuclotide sequence, 
        iterations: number of simulations to run
    
    Returns:
        a dictionary containing P values, and distances for missense and nonsense,
        de novos. Missing data is represented by "NA".
    """
    
    if len(gene.transcripts) == 0:
        nan = float('nan')
        return {'miss_dist': nan, 'miss_prob': nan, 'nons_prob': nan, 'nons_dist': nan}
    
    transcript = gene.canonical
    
    missense_events = get_de_novos_in_transcript(transcript, de_novos["missense"])
    nonsense_events = get_de_novos_in_transcript(transcript, de_novos["nonsense"])
    
    # The rates either are passed in as lists of lists of sequence-context 
    # based rates which can be applied to a given sequence, or as 
    # dictionaries of rates per position/allele from genome based reference 
    # sources. If the latter, then we have to load the rates within the 
    # transcript region
    _rates = mut_dict
    if not isinstance(mut_dict, list):
        _rates = mut_dict(transcript.chrom, transcript.start, transcript.end)
    
    rates = SiteRates(transcript, _rates)
    
    (miss_dist, miss_prob) = get_structure_p_value(transcript, rates, structure, iterations, "missense", missense_events)
    (nons_dist, nons_prob) = get_structure_p_value(transcript, rates, structure, iterations, "lof", nonsense_events)
    
    return {"miss_dist": miss_dist,
            "nons_dist": nons_dist,
            "miss_prob": miss_prob,
            "nons_prob": nons_prob,
    }
    
