""" Script to investigate the probability of multiple mutations clustering
within a single gene.
"""

from __future__ import print_function, division, absolute_import

import sys
import os
import argparse
import math

import scipy.stats

from denovonear.load_gene import get_deprecated_gene_ids, load_gene, \
    get_de_novos_in_transcript
from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_mutation_rates import load_trincleotide_mutation_rates
from denovonear.load_de_novos import load_de_novos
from denovonear.site_specific_rates import SiteRates
from denovonear.analyse_de_novo_clustering import AnalyseDeNovoClustering

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
    parser.add_argument("--rates", dest="rates_path", required=True, \
        help="Path to file containing trinucleotide mutation rates.")
    parser.add_argument("--deprecated-genes", dest="deprecated_genes_path", \
        help="deprecated gene IDs filename")
    parser.add_argument("--genome-build", dest="genome_build", choices=["grch37",
        "GRCh37", "grch38", "GRCh38"], default="grch37", help="Genome build "+ \
        "that the de novo coordinates are based on (GrCh37 or GRCh38")
    parser.add_argument("--cache-folder", dest="cache_dir", \
        default=os.path.join(os.path.dirname(__file__), "cache"), help="folder \
        to cache Ensembl data into (defaults to clustering code directory)")
    
    args = parser.parse_args()
    
    args.genome_build = args.genome_build.lower()
    
    return args

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

def analyse_gene(gene_id, iterations, ensembl, de_novos, old_gene_ids, mut_dict):
    """ run the analysis code for a single gene
    
    Args:
        gene_id: HGNC symbol for a gene
        iterations: number of simulations to run
        ensembl: EnsemblRequest object, for obtaing info from ensembl
        known_de_novos: dictionary of de novo positions for the HGNC gene,
            indexed by functional type
        old_gene_ids: dictionary of updated HGNC symbols, indexed by their old ID
        mut_dict: dictionary of mutation rates, indexed by trinuclotide sequence
    
    Returns:
        a dictionary containing P values, and distances for missense, nonsense,
        and synonymous de novos events. Missing data is represented by "NA".
    """
    
    # fix HGNC IDs that have been discontinued in favour of other gene IDs
    if gene_id in old_gene_ids:
        gene_id = old_gene_ids[gene_id]
    
    missense = de_novos["missense"]
    nonsense = de_novos["nonsense"]
    
    # load the set of transcripts that are the  minimum set of transcripts
    # required to contain all the de novos, unless we can't find any coding
    # transcripts that contain the de novos.
    try:
        transcripts = load_gene(ensembl, gene_id, missense + nonsense)
    except IndexError as e:
        print(e)
        return None
    
    probs = {"miss_prob": [], "nons_prob": []}
    dists = {"miss_dist": [], "nons_dist": []}
    
    for transcript in transcripts:
        
        missense_events = get_de_novos_in_transcript(transcript, missense)
        nonsense_events = get_de_novos_in_transcript(transcript, nonsense)
        
        site_weights = SiteRates(transcript, mut_dict)
        
        print("simulating clustering")
        clust = AnalyseDeNovoClustering(transcript, site_weights, iterations)
        
        (miss_dist, miss_prob) = clust.analyse_missense_and_splice_region(missense_events)
        (nons_dist, nons_prob) = clust.analyse_lof(nonsense_events)
        
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
        dists[key] = ",".join(dists[key])
    
    probs = combine_p_values(probs)
    probs.update(dists)
    
    return probs

def main():
    
    args = get_options()
    
    # load all the data
    ensembl = EnsemblRequest(args.cache_dir, args.genome_build)
    mut_dict = load_trincleotide_mutation_rates(args.rates_path)
    
    old_gene_ids = {}
    # only load the old gene ID converter if we have specified the file
    if args.deprecated_genes_path is not None:
        old_gene_ids = get_deprecated_gene_ids(args.deprecated_genes_path)
    
    known_de_novos = load_de_novos(args.input)
    
    output = open(args.output, "w")
    output.write("\t".join(["gene_id", "mutation_category", "events_n", \
        "dist", "probability"]) + "\n")
    
    initial_iterations = 1000000
    for gene_id in sorted(known_de_novos):
        
        de_novos = known_de_novos[gene_id]
        
        if len(de_novos["missense"] + de_novos["nonsense"]) < 2:
            continue
        
        probs = analyse_gene(gene_id, initial_iterations, ensembl, de_novos, \
            old_gene_ids, mut_dict)
        
        if probs is None:
            continue
        
        output.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(gene_id, "missense", \
            len(de_novos["missense"]), probs["miss_dist"], probs["miss_prob"]))
        output.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(gene_id, "nonsense", \
            len(de_novos["nonsense"]), probs["nons_dist"], probs["nons_prob"]))
        
        # sys.exit()

if __name__ == '__main__':
    main()
