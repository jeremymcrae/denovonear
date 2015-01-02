""" Script to investigate the probability of multiple mutations clustering 
within a single gene.
"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import os
import argparse

from src.load_gene import get_deprecated_gene_ids, load_gene, load_conservation
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
        iterations = initial_iterations
        # gene_id = "PACS1"
        # print(gene_id)
        
        func_events = known_de_novos[gene_id]["functional"]
        missense_events = known_de_novos[gene_id]["missense"]
        nonsense_events = known_de_novos[gene_id]["nonsense"]
        synonymous_events = known_de_novos[gene_id]["synonymous"]
        
        # don't analyse genes with only one de novo functional mutation, and 
        # for now, exclude genes with numerous events
        if len(synonymous_events) < 2 and len(missense_events) < 2 and len(nonsense_events) < 2:
            continue
        
        # fix HGNC IDs that have been discontinued in favour of other gene IDs
        if gene_id in old_gene_ids:
            gene_id = old_gene_ids[gene_id]
        
        try:
            transcript = load_gene(ensembl, gene_id, func_events)
            # transcript = load_conservation(transcript, conservation_folder)
        except IndexError as e:
            print(e)
            continue
        
        site_weights = SiteRates(transcript, mut_dict, use_coverage=use_coverage)
        if coverage_dir is not None:
            site_weights.set_coverage_dir(coverage_dir)
        
        print("simulating clustering")
        probs = AnalyseDeNovoClustering(transcript, site_weights, iterations)
        
        (miss_dist, miss_prob) = probs.analyse_missense_and_splice_region(missense_events)
        (nons_dist, nons_prob) = probs.analyse_lof(nonsense_events)
        (syn_dist, syn_prob) = probs.analyse_synonymous(synonymous_events)
        
        # [cons_func, cons_func_p, cons_miss, cons_miss_p, cons_nons, \
        #     cons_nons_p] = ["NA"] * 6
        # if hasattr(transcript, "conservation_scores"):
        #     probs = AnalyseDeNovoConservation(transcript, site_weights, iterations)
            
        #     # (cons_func, cons_func_p) = probs.analyse_functional(func_events)
        #     (cons_miss, cons_miss_p) = probs.analyse_missense(missense_events)
        #     (cons_nons, cons_nons_p) = probs.analyse_nonsense(nonsense_events)
        
        output.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(gene_id, "missense", \
            len(missense_events), miss_dist, miss_prob))
        output.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(gene_id, "nonsense", \
            len(nonsense_events), nons_dist, nons_prob))
        output.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(gene_id, "synonymous", \
            len(synonymous_events), syn_dist, syn_prob))
        
        # sys.exit()

if __name__ == '__main__':
    main()





