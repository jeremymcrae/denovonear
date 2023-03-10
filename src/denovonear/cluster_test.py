
from math import log, isnan

from scipy.stats import chi2

from denovonear.load_gene import load_gene, get_de_novos_in_transcript, \
    minimise_transcripts
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.load_de_novos import load_de_novos
from denovonear.site_specific_rates import SiteRates
from denovonear.simulate import get_p_value

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

def cluster_de_novos(symbol, de_novos, gene, iterations=1000000, 
        mut_dict=None):
    """ analysis proximity cluster of de novos in a single gene
    
    Args:
        symbol: HGNC symbol for a gene
        de_novos: dictionary of de novo positions for the HGNC gene,
        indexed by functional type
        gene: denovnonear.gencode.Gene object
        iterations: number of simulations to run
        mut_dict: dictionary of mutation rates, indexed by trinuclotide sequence
    
    Returns:
        a dictionary containing P values, and distances for missense, nonsense,
        and synonymous de novos events. Missing data is represented by "NA".
    """
    
    if mut_dict is None:
        mut_dict = load_mutation_rates()
    
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
    transcripts = [x for x in transcripts if x.get_name() in minimized]
    
    if len(transcripts) == 0:
        nan = float('nan')
        return {'miss_dist': nan, 'miss_prob': nan, 'nons_prob': nan, 'nons_dist': nan}
    
    probs = {"miss_prob": [], "nons_prob": []}
    dists = {"miss_dist": [], "nons_dist": []}
    
    for transcript in transcripts:
        
        missense_events = get_de_novos_in_transcript(transcript, missense)
        nonsense_events = get_de_novos_in_transcript(transcript, nonsense)
        
        rates = SiteRates(transcript, mut_dict)
        
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
