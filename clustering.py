""" Script to investigate the probability of multiple mutations clustering 
within a single gene.
"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import os
import argparse

from src.load_gene import get_deprecated_gene_ids, get_transcript_lengths, \
    construct_gene_object, check_denovos_in_gene, load_gene, load_conservation
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
    
    parser = argparse.ArgumentParser(description='examine mutation clustering in genes')
    parser.add_argument("--in", dest="input", help="input filename for file listing known mutations in genes")
    parser.add_argument("--out", dest="output", help="output filename")
    parser.add_argument("--rates", dest="mut_rates", help="mutation rates filename")
    parser.add_argument("--deprecated-genes", dest="deprecated_genes", help="deprecated gene IDs filename")
    
    args = parser.parse_args()
    
    return args.input, args.output, args.mut_rates, args.deprecated_genes


def main():
    
    input_file, output_file, mut_rates_file, deprecated_gene_id_file = get_options()
    
    # load all the data
    ensembl = EnsemblRequest()
    mut_dict = load_trincleotide_mutation_rates(mut_rates_file)
    old_gene_ids = get_deprecated_gene_ids(deprecated_gene_id_file)
    known_de_novos = load_known_de_novos(input_file)
    
    output = open(output_file, "w")
    output.write("\t".join(["gene_id", \
        "functional_events_n", "functional_dist", "functional_probability", 
        "missense_events_n", "missense_dist", "missense_probability", 
        "nonsense_events_n", "nonsense_distance", "nonsense_dist_probability"]) + "\n")
    
    initial_iterations = 1000000
    for gene_id in known_de_novos:
        iterations = initial_iterations
        print(gene_id)
        
        func_events = known_de_novos[gene_id]["functional"]
        missense_events = known_de_novos[gene_id]["missense"]
        nonsense_events = known_de_novos[gene_id]["nonsense"]
        
        # don't analyse genes with only one de novo functional mutation
        if len(func_events) < 2:
            continue
        
        # fix HGNC IDs that have been discontinued in favour of other gene IDs
        if gene_id in old_gene_ids:
            gene_id = old_gene_ids[gene_id]
        
        try:
            transcript = load_gene(ensembl, gene_id, func_events)
            # transcript = load_conservation(transcript, conservation_folder)
        except IndexError:
            continue
        
        site_weights = SiteRates(transcript, mut_dict)
        
        print("simulating clustering")
        probs = AnalyseDeNovoClustering(transcript, site_weights, iterations)
        
        (func_dist, func_prob) = probs.analyse_functional(func_events)
        (miss_dist, miss_prob) = probs.analyse_missense(missense_events)
        (nons_dist, nons_prob) = probs.analyse_nonsense(nonsense_events)
        
        # [cons_func, cons_func_p, cons_miss, cons_miss_p, cons_nons, \
        #     cons_nons_p] = ["NA"] * 6
        # if hasattr(transcript, "conservation_scores"):
        #     probs = AnalyseDeNovoConservation(transcript, site_weights, iterations)
            
        #     # (cons_func, cons_func_p) = probs.analyse_functional(func_events)
        #     (cons_miss, cons_miss_p) = probs.analyse_missense(missense_events)
        #     (cons_nons, cons_nons_p) = probs.analyse_nonsense(nonsense_events)
        
        # output.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n".\
        #     format(gene_id, \
        #     len(missense_events), miss_dist, miss_prob, cons_miss, cons_miss_p, \
        #     len(nonsense_events), nons_dist, nons_prob, cons_nons, cons_nons_p ))
        
        output.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(gene_id, \
            len(func_events), func_dist, func_prob, \
            len(missense_events), miss_dist, miss_prob, \
            len(nonsense_events), nons_dist, nons_prob ))
        
        # sys.exit()

if __name__ == '__main__':
    main()





