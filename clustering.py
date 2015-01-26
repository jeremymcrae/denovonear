""" Script to investigate the probability of multiple mutations clustering 
within a single gene.
"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import os
import argparse
import math

import scipy.stats

from src.load_gene import get_deprecated_gene_ids, load_gene, load_conservation, \
    get_de_novos_in_transcript
from src.ensembl_requester import EnsemblRequest
from src.load_mutation_rates import load_trincleotide_mutation_rates
from src.load_known_de_novos import load_known_de_novos
from src.load_conservation_scores import load_conservation_scores
from src.site_specific_rates import SiteRates
from src.analyse_de_novo_clustering import AnalyseDeNovoClustering
from src.analyse_de_novo_conservation import AnalyseDeNovoConservation

def get_options():
    """ get the command line switches
    """
    
    parser = argparse.ArgumentParser(description="Examines the proximity \
        clustering of de novo mutations in genes.")
    parser.add_argument("--in", dest="input", required=True, help="Path to \
        file listing known mutations in genes. See example file in data folder \
        for format.")
    parser.add_argument("--out", dest="output", required=True, help="output \
        filename")
    parser.add_argument("--rates", dest="mut_rates", required=True, \
        help="Path to file containing trinucleotide mutation rates.")
    parser.add_argument("--deprecated-genes", dest="deprecated_genes", \
        help="deprecated gene IDs filename")
    parser.add_argument("--genome-build", dest="genome_build", choices=["grch37",
        "GRCh37", "grch38", "GRCh38"], default="grch37", help="Genome build "+ \
        "that the de novo coordinates are based on (GrCh37 or GRCh38")
    parser.add_argument("--cache-folder", dest="cache_folder", \
        default=os.path.join(os.path.dirname(__file__), "cache"), help="folder \
        to cache Ensembl data into (defaults to clustering code directory)")
    parser.add_argument("--coverage-adjust", default=False, action="store_true", \
        help="whether to adjust site mutation rates for sequencing coverage.")
    parser.add_argument("--coverage-dir", help="location of ExAC coverage files")
    
    args = parser.parse_args()
    
    return args.input, args.output, args.mut_rates, args.deprecated_genes, \
        args.cache_folder, args.genome_build.lower(), args.coverage_adjust, \
        args.coverage_dir

def combine_p_values(probs):
    """ Combine the P values from different transcripts.
    
    This returns P values for each mutation type for a gene. We occasionally
    have multiple P values for a mutation type, obtained from different
    transcripts for the gene. If we have only one P value for the gene for the
    mutation type, we just use that value, if we don't have any data, we use
    "NA", otherwise we combine the P values from different transcripts using
    Fisher's combined test.
    
    Args:
        probs: Dictionary of lists of P values from different transcripts,
            indexed by functional type.
    
    Returns:
        Dictionary of P values, indexed by the mutation type
    """
    
    fixed_probs = {}
    for key in probs:
        values = probs[key]
        
        # drop out the NA values
        values = [x for x in values if x != "NA"]
        
        if len(values) == 0:
            fixed_probs[key] = "NA"
        elif len(values) == 1:
            fixed_probs[key] = values[0]
        else:
            # use Fisher's combined method to estimate the P value from multiple
            # P values. The chi square statistic is -2*sum(ln(P values))
            values = [math.log(x) for x in values]
            chi_square = -2 * sum(values)
            df = 2 * len(values)
            
            # estimate the P value using the chi square statistic and degrees of
            # freedom
            p_value = 1 - scipy.stats.chi2.cdf(chi_square, df)
            fixed_probs[key] = p_value
    
    return fixed_probs

def analyse_gene(gene_id, iterations, ensembl, de_novos, old_gene_ids, mut_dict, use_coverage, coverage_dir):
    """ run the analysis code for a single gene
    
    Args:
        gene_id: HGNC symbol for a gene
        iterations: number of simulations to run
        ensembl: EnsemblRequest object, for obtaing info from ensembl
        known_de_novos: dictionary of de novo positions for the HGNC gene, 
            indexed by functional type
        old_gene_ids: dictionary of updated HGNC symbols, indexed by their old ID
        mut_dict: dictionary of mutation rates, indexed by trinuclotide sequence
        use_coverage: whether to compensate for the depth of sequence coverage
        coverage_dir: path to folder containing coverage data, or None
    
    Returns:
        a dictionary containing P values, and distances for missense, nonsense, 
        and synonymous de novos events. Missing data is represented by "NA".
    """
    
    # fix HGNC IDs that have been discontinued in favour of other gene IDs
    if gene_id in old_gene_ids:
        gene_id = old_gene_ids[gene_id]
    
    func_events = de_novos["functional"]
    missense = de_novos["missense"]
    nonsense = de_novos["nonsense"]
    synonymous = de_novos["synonymous"]
    
    # load the set of transcripts that are the  minimum set of transcripts 
    # required to contain all the de novos, unless we can't find any coding
    # transcripts that contain the de novos.
    try:
        transcripts = load_gene(ensembl, gene_id, func_events)
    except IndexError as e:
        print(e)
        return None
    
    probs = {"miss_prob": [], "nons_prob": [], "syn_prob": []}
    dists = {"miss_dist": [], "nons_dist": [], "syn_dist": []}
    
    for transcript in transcripts:
        
        missense_events = get_de_novos_in_transcript(transcript, missense)
        nonsense_events = get_de_novos_in_transcript(transcript, nonsense)
        synonymous_events = get_de_novos_in_transcript(transcript, synonymous)
        
        site_weights = SiteRates(transcript, mut_dict, use_coverage=use_coverage)
        if coverage_dir is not None:
            site_weights.set_coverage_dir(coverage_dir)
        
        print("simulating clustering")
        clust = AnalyseDeNovoClustering(transcript, site_weights, iterations)
    
        (miss_dist, miss_prob) = clust.analyse_missense_and_splice_region(missense_events)
        (nons_dist, nons_prob) = clust.analyse_lof(nonsense_events)
        (syn_dist, syn_prob) = clust.analyse_synonymous(synonymous_events)
        
        dists["miss_dist"].append(miss_dist)
        dists["nons_dist"].append(nons_dist)
        dists["syn_dist"].append(syn_dist)
        probs["miss_prob"].append(miss_prob)
        probs["nons_prob"].append(nons_prob)
        probs["syn_prob"].append(syn_prob)
        
        # remove the de novos analysed in the current transcript, so that 
        # analysis of subsequent transcripts uses independent events. NOTE THAT
        # THIS MIGHT MISS SOME CLUSTERING ACROSS MUTUALLY EXCLUSIVE TRANSCRIPTS
        # IF THE DE NOVO EVENTS ARE NEAR THE TRANSCRIPT DIVERGENCE.
        missense = [x for x in missense if x not in missense_events]
        nonsense = [x for x in nonsense if x not in  nonsense_events]
        synonymous = [x for x in synonymous if x not in  synonymous_events]
        
    for key in dists:
        dists[key] = ",".join(dists[key])
    
    probs = combine_p_values(probs)
    probs.update(dists)
    
    return probs  

def main():
    
    input_file, output_file, rates_file, old_gene_id_file, cache_dir, \
        genome_build, use_coverage, coverage_dir = get_options()
    
    # load all the data
    ensembl = EnsemblRequest(cache_dir, genome_build)
    mut_dict = load_trincleotide_mutation_rates(rates_file)
    
    old_gene_ids = {}
    # only load the old gene ID converter if we have specified the file
    if old_gene_id_file is not None:
        old_gene_ids = get_deprecated_gene_ids(old_gene_id_file)
    
    known_de_novos = load_known_de_novos(input_file)
    
    output = open(output_file, "w")
    output.write("\t".join(["gene_id", "mutation_category", "events_n", \
        "dist", "probability"]) + "\n")
    
    initial_iterations = 1000000
    for gene_id in sorted(known_de_novos):
        
        de_novos = known_de_novos[gene_id]
        probs = analyse_gene(gene_id, initial_iterations, ensembl, de_novos, old_gene_ids, mut_dict, use_coverage, coverage_dir)
        
        if probs is None:
            continue
        
        output.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(gene_id, "missense", \
            len(de_novos["missense"]), probs["miss_dist"], probs["miss_prob"]))
        output.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(gene_id, "nonsense", \
            len(de_novos["nonsense"]), probs["nons_dist"], probs["nons_prob"]))
        output.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(gene_id, "synonymous", \
            len(de_novos["synonymous"]), probs["syn_dist"], probs["syn_prob"]))
        
        # sys.exit()

if __name__ == '__main__':
    main()





