""" Script to identify transcripts containing de novos
"""

import os
import argparse

from denovonear.load_gene import count_de_novos_per_transcript, minimise_transcripts
from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_de_novos import load_de_novos

def get_options():
    """ get the command line switches
    """
    
    parser = argparse.ArgumentParser(description="Identify transcripts "
        "for a gene containing de novo events.")
    parser.add_argument("--de-novos", required=True, help="Path to "
        "file listing de novo variants in genes.")
    parser.add_argument("--out", required=True, help="output filename")
    parser.add_argument("--genome-build", choices=["grch37", "GRCh37", "grch38",
        "GRCh38"], default="grch37", help="Genome build that the de novo "
            "coordinates are based on (GrCh37 or GRCh38")
    parser.add_argument("--cache-folder",
        default=os.path.join(os.path.dirname(__file__), "cache"), help="folder "
            "to cache Ensembl data into (defaults to clustering code directory)")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--all-transcripts", action="store_true", default=False,
        help="Flag if you want to identify all transcripts with more than "
        "one de novo on it.")
    group.add_argument("--minimise-transcripts", action="store_true",
        default=False, help="Flag if you want to identify the minimal set of "
        "transcripts to contain all de novos.")
    
    return parser.parse_args()

def main():
    
    args = get_options()
    ensembl = EnsemblRequest(args.cache_folder, args.genome_build)
    de_novos = load_de_novos(args.de_novos, exclude_indels=False)
    
    output = open(args.out, "w")
    output.write("hgnc_symbol\ttranscript_id\tlength\tde_novos\n")
    
    for symbol in sorted(de_novos):
        print(symbol)
        func_events = de_novos[symbol]["missense"] + de_novos[symbol]["nonsense"]
        
        # find the counts per transcript, depending on whether we want to count
        # for all transcripts containing one or more de novos, or to find the
        # minimum set of transcripts to contain the de novos
        try:
            if args.all_transcripts:
                counts = count_de_novos_per_transcript(ensembl, symbol, func_events)
            elif args.minimal_transcripts:
                counts = minimise_transcripts(ensembl, symbol, func_events)
        except (ValueError, IndexError):
            print("error occured with {0}".format(symbol))
            continue
        
        # write the transcript details to a file
        for key in counts:
            line = "{}\t{}\t{}\t{}\n".format(symbol, key, counts[key]["len"],
                counts[key]["n"])
            output.write(line)
        
    output.close()

if __name__ == '__main__':
    main()
