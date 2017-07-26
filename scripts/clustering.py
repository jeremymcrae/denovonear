""" Script to investigate the probability of multiple mutations clustering
within a single gene.
"""

from __future__ import print_function, division, absolute_import

import os
import argparse

from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.load_de_novos import load_de_novos
from denovonear.cluster_test import cluster_de_novos

def get_options():
    """ get the command line switches
    """
    
    parser = argparse.ArgumentParser(description="Examines the proximity" \
        "clustering of de novo mutations in genes.")
    parser.add_argument("--in", dest="input", required=True, help="Path to" \
        "file listing known mutations in genes. See example file in data folder" \
        "for format.")
    parser.add_argument("--out", dest="output", required=True, help="output \
        filename")
    parser.add_argument("--rates",
        help="Path to file containing sequence context-based mutation rates.")
    parser.add_argument("--deprecated-genes", dest="deprecated_genes_path",
        help="deprecated gene IDs filename")
    parser.add_argument("--genome-build", dest="genome_build", choices=["grch37",
        "GRCh37", "grch38", "GRCh38"], default="grch37", help="Genome build " \
        "that the de novo coordinates are based on (GrCh37 or GRCh38")
    parser.add_argument("--cache-folder", dest="cache_dir",
        default=os.path.join(os.path.dirname(__file__), "cache"), help="folder" \
        "to cache Ensembl data into (defaults to clustering code directory)")
    
    args = parser.parse_args()
    
    args.genome_build = args.genome_build.lower()
    
    return args

def main():
    
    args = get_options()
    
    # load all the data
    ensembl = EnsemblRequest(args.cache_dir, args.genome_build)
    mut_dict = load_mutation_rates(args.rates)
    
    old_gene_ids = {}
    # only load the old gene ID converter if we have specified the file
    if args.deprecated_genes_path is not None:
        old_gene_ids = get_deprecated_gene_ids(args.deprecated_genes_path)
    
    known_de_novos = load_de_novos(args.input)
    
    output = open(args.output, "w")
    output.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")
    
    iterations = 1000000
    for symbol in sorted(known_de_novos):
        
        de_novos = known_de_novos[symbol]
        
        if len(de_novos["missense"] + de_novos["nonsense"]) < 2:
            continue
        
        # fix HGNC IDs that have been discontinued in favour of other gene IDs
        if symbol in old_gene_ids:
            symbol = old_gene_ids[gene_id]
        
        probs = cluster_de_novos(symbol, de_novos, iterations, ensembl, mut_dict)
        
        if probs is None:
            continue
        
        output.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(symbol, "missense", \
            len(de_novos["missense"]), probs["miss_dist"], probs["miss_prob"]))
        output.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(symbol, "nonsense", \
            len(de_novos["nonsense"]), probs["nons_dist"], probs["nons_prob"]))

if __name__ == '__main__':
    main()
