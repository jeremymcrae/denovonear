""" Script to investigate the probability of multiple mutations clustering 
within a single gene.
"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import os
import math
import argparse

from src.load_gene import construct_gene_object
from src.ensembl_requester import EnsemblRequest
from src.load_mutation_rates import load_trincleotide_mutation_rates
from src.site_specific_rates import SiteRates

def get_options():
    """ get the command line switches
    """
    
    parser = argparse.ArgumentParser(description="determine mutation rates \
        for genes given transcript IDs.")
    parser.add_argument("--transcripts", dest="input", required=True, help="Path to \
        file listing transcript IDs.")
    parser.add_argument("--out", dest="output", required=True, help="output \
        filename")
    parser.add_argument("--rates", dest="mut_rates", required=True, \
        help="Path to file containing trinucleotide mutation rates.")
    parser.add_argument("--deprecated-genes", dest="deprecated_genes", \
        help="deprecated gene IDs filename")
    parser.add_argument("--cache-folder", dest="cache_folder", \
        default=os.path.join(os.path.dirname(__file__), "cache"), help="folder \
        to cache Ensembl data into (defaults to clustering code directory)")
    
    args = parser.parse_args()
    
    return args.input, args.output, args.mut_rates, args.deprecated_genes, \
        args.cache_folder

def load_transcripts(path):
    """ load a file listing transcript IDs per line
    
    Args:
        path: path to file containing transcript IDs, one per line
    
    Returns:
        list of transcript IDs eg ["ENST00000315684", "ENST00000485511"]
    """
    
    transcript_ids = []
    with open(path, "r") as f:
        for line in f:
            transcript_ids.append(line.strip())
            
    return transcript_ids

def main():
    
    input_file, output_file, rates_file, old_gene_id_file, cache_dir = get_options()
    
    # load all the data
    ensembl = EnsemblRequest(cache_dir)
    mut_dict = load_trincleotide_mutation_rates(rates_file)
    
    transcript_ids = load_transcripts(input_file)
    
    output = open(output_file, "w")
    output.write("transcript_id\tmissense_rate\tnonsense_rate\n")
    
    for transcript_id in transcript_ids:
        
        transcript = construct_gene_object(ensembl, transcript_id)
        site_weights = SiteRates(transcript, mut_dict)
        
        missense = math.log10(site_weights.get_missense_rates_for_gene().cum_probs[-1])
        nonsense = math.log10(site_weights.get_nonsense_rates_for_gene().cum_probs[-1])
        
        line = "{0}\t{1}\t{2}\n".format(transcript_id, missense, nonsense)
        output.write(line)
        
    output.close()

if __name__ == '__main__':
    main()





