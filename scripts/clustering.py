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
    
    parser = argparse.ArgumentParser(description="Examines the proximity "
        "clustering of de novo mutations in genes.")
    parser.add_argument("--in", dest="input", required=True, help="Path to "
        "file listing known mutations in genes. See example file in data folder "
        "for format.")
    parser.add_argument("--out", required=True, help="output filename")
    parser.add_argument("--rates",
        help="Path to file containing sequence context-based mutation rates.")
    parser.add_argument("--genome-build", choices=["grch37",
        "GRCh37", "grch38", "GRCh38"], default="grch37", help="Genome build "
        "that the de novo coordinates are based on (GrCh37 or GRCh38")
    parser.add_argument("--cache-folder",
        default=os.path.join(os.path.dirname(__file__), "cache"), help="folder "
        "to cache Ensembl data into (defaults to clustering code directory)")
    
    return parser.parse_args()

def main():
    
    args = get_options()
    
    # load all the data
    ensembl = EnsemblRequest(args.cache_folder, args.genome_build.lower())
    mut_dict = load_mutation_rates(args.rates)
    de_novos = load_de_novos(args.input)
    
    output = open(args.output, "w")
    output.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")
    
    iterations = 1000000
    for symbol in sorted(de_novos):
        
        de_novos = de_novos[symbol]
        
        if len(de_novos["missense"] + de_novos["nonsense"]) < 2:
            continue
        
        probs = cluster_de_novos(symbol, de_novos, iterations, ensembl, mut_dict)
        
        if probs is None:
            continue
        
        output.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense",
            len(de_novos["missense"]), probs["miss_dist"], probs["miss_prob"]))
        output.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "nonsense",
            len(de_novos["nonsense"]), probs["nons_dist"], probs["nons_prob"]))

if __name__ == '__main__':
    main()
