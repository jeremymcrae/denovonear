""" Script to identify transcripts containing de novos
"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import os
import math
import argparse

from src.load_gene import construct_gene_object, get_deprecated_gene_ids, \
    count_de_novos_per_transcript, minimise_transcripts
from src.ensembl_requester import EnsemblRequest
from src.load_mutation_rates import load_trincleotide_mutation_rates
from src.site_specific_rates import SiteRates
from src.load_known_de_novos import load_known_de_novos

def get_options():
    """ get the command line switches
    """
    
    parser = argparse.ArgumentParser(description="Identify transcripts \
        for a gene containing de novo events.")
    parser.add_argument("--de-novos", dest="input", required=True, help="Path to \
        file listing de novo variants in genes.")
    parser.add_argument("--out", dest="output", required=True, help="output \
        filename")
    parser.add_argument("--deprecated-genes", dest="deprecated_genes", \
        help="deprecated gene IDs filename")
    parser.add_argument("--genome-build", dest="genome_build", choices=["grch37",
        "GRCh37", "grch38", "GRCh38"], default="grch37", help="Genome build "+ \
        "that the de novo coordinates are based on (GrCh37 or GRCh38")
    parser.add_argument("--cache-folder", dest="cache_folder", \
        default=os.path.join(os.path.dirname(__file__), "cache"), help="folder \
        to cache Ensembl data into (defaults to clustering code directory)")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--all-transcripts", action="store_true", default=False,\
        help="Flag if you want to identify all transcripts with more than" + \
        "one de novo on it.")
    group.add_argument("--minimise-transcripts", action="store_true", \
        default=False, help="Flag if you want to identify the minimal set of" +\
        "transcripts to contain all de novos.")
    
    args = parser.parse_args()
    
    return args.input, args.output, args.deprecated_genes, \
        args.cache_folder, args.genome_build.lower(), args.all_transcripts, \
        args.minimise_transcripts

def main():
    
    input_file, output_file, old_gene_id_file, cache_dir, genome_build, \
        all_transcripts, minimal_transcripts = get_options()
    
    # load all the data
    ensembl = EnsemblRequest(cache_dir, genome_build)
    
    old_gene_ids = {}
    if old_gene_id_file is not None:
        old_gene_ids = get_deprecated_gene_ids(old_gene_id_file)
    
    known_de_novos = load_known_de_novos(input_file)
    
    output = open(output_file, "w")
    output.write("hgnc_symbol\ttranscript_id\tlength\tde_novos\n")
    
    for gene_id in sorted(known_de_novos):
        func_events = known_de_novos[gene_id]["functional"]
        
        # fix HGNC IDs that have been discontinued in favour of other gene IDs
        if gene_id in old_gene_ids:
            gene_id = old_gene_ids[gene_id]
        
        # find the counts per transcript, depending on whether we want to count
        # for all transcripts containing one or more de novos, or to find the 
        # minimum set of transcripts to contain the de novos
        if all_transcripts:
            counts = count_de_novos_per_transcript(ensembl, gene_id, func_events)
        elif minimal_transcripts:
            counts = minimise_transcripts(ensembl, gene_id, func_events)
        
        # write the transcript details to a file
        for (transcript_id, count, length) in counts:
            line = "{0}\t{1}\t{2}\t{3}\n".format(gene_id, transcript_id, length, count)
            output.write(line)
        
    output.close()

if __name__ == '__main__':
    main()


