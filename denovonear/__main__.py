""" Script to investigate the probability of multiple mutations clustering
within a single gene.
"""

import os
import sys
import argparse

from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.load_de_novos import load_de_novos
from denovonear.cluster_test import cluster_de_novos

from denovonear.load_gene import (construct_gene_object,
    count_de_novos_per_transcript, minimise_transcripts)
from denovonear.site_specific_rates import SiteRates
from denovonear.frameshift_rate import include_frameshift_rates
from denovonear.log_transform_rates import log_transform

def clustering(ensembl, mut_dict, output, args):
    
    de_novos = load_de_novos(args.input)
    
    output.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")
    
    iterations = 1000000
    for symbol in sorted(de_novos):
        
        if len(de_novos[symbol]["missense"] + de_novos[symbol]["nonsense"]) < 2:
            continue
        
        probs = cluster_de_novos(symbol, de_novos[symbol], iterations, ensembl, mut_dict)
        
        if probs is None:
            continue
        
        output.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense",
            len(de_novos[symbol]["missense"]), probs["miss_dist"], probs["miss_prob"]))
        output.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "nonsense",
            len(de_novos[symbol]["nonsense"]), probs["nons_dist"], probs["nons_prob"]))

def find_transcripts(ensembl, mut_dict, output, args):
    
    de_novos = load_de_novos(args.de_novos)
    
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

def load_genes(path):
    """ load a file listing gene and transcript IDs
    
    Args:
        path: path to file containing gene IDs and transcript IDs e.g.
            gene_1    transcript_1.1    length_1    denovo_count
            gene_2    transcript_2.1    length_3    denovo_count
    
    Returns:
        dict of transcripts eg {'CTC1': ["ENST00000315684", "ENST00000485511"]}
    """
    
    with open(path, 'rt') as f:
        lines = [ x.split('\t')[:2] for x in f if not x.startswith('hgnc') ]
    
    transcripts = {}
    for symbol, tx in lines:
        if symbol not in transcripts:
            transcripts[symbol] = []
        transcripts[symbol].append(tx)
    
    return transcripts

def get_mutation_rates(transcripts, mut_dict, ensembl):
    """ determines mutation rates per functional category for transcripts
    
    Args:
        transcripts: list of transcript IDs for a gene
        mut_dict: dictionary of local sequence context mutation rates
        ensembl: EnsemblRequest object, to retrieve information from Ensembl.
    
    Returns:
        tuple of (rates, merged transcript, and transcript CDS length)
    """
    
    rates = {'missense': 0, 'nonsense': 0, 'splice_lof': 0,
        'splice_region': 0, 'synonymous': 0}
    combined = None
    
    for tx_id in transcripts:
        try:
            tx = construct_gene_object(ensembl, tx_id)
        except ValueError:
            continue
        
        if len(tx.get_cds_sequence()) % 3 != 0:
            raise ValueError("anomalous_coding_sequence")
        
        # ignore mitochondrial genes
        if tx.get_chrom() == "MT":
            continue
        
        sites = SiteRates(tx, mut_dict, masked_sites=combined)
        combined = tx + combined
        
        for cq in ['missense', 'nonsense', 'splice_lof', 'splice_region', 'synonymous']:
            rates[cq] += sites[cq].get_summed_rate()
    
    if combined is None:
        raise ValueError('no tx found')
    
    length = combined.get_coding_distance(combined.get_cds_start(),
        combined.get_cds_end())
    
    return rates, combined, length

def gene_rates(ensembl, mut_dict, output, args):
    
    transcripts = load_genes(args.genes)
    
    header = ['transcript_id', 'chrom', 'length', 'missense_rate', 'nonsense_rate',
        'splice_lof_rate', 'splice_region_rate', 'synonymous_rate']
    output.write('\t'.join(header) + "\n")
    
    for symbol in sorted(transcripts):
        print(symbol)
        try:
            rates, tx, length = get_mutation_rates(transcripts[symbol],
                mut_dict, ensembl)
            # log transform rates, for consistency with Samocha et al.
            line = "{}\t{}\t{}\t{}".format(symbol, tx.get_chrom(), length, log_transform(rates))
        except (ValueError, KeyError) as error:
            print("{}\t{}\n".format(symbol, error))
            line = "{}\t{}\tNA\tNA\tNA\tNA\tNA\tNA".format(symbol, tx.get_chrom())
        
        output.write(line + '\n')
    
    output.close()
    include_frameshift_rates(args.out)

def get_options():
    """ get the command line switches
    """
    
    parser = argparse.ArgumentParser(description='denovonear cli interface')
    
    ############################################################################
    # CLI options in common
    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument("--out", help="output filename")
    parent.add_argument("--rates",
        help="Path to file containing sequence context-based mutation rates.")
    parent.add_argument("--genome-build", choices=["grch37",
        "GRCh37", "grch38", "GRCh38"], default="grch37", help="Genome build "
        "that the de novo coordinates are based on (GRCh37 or GRCh38")
    parent.add_argument("--cache-folder",
        default=os.path.join(os.path.expanduser('~'), ".cache", 'denovonear'),
        help="where to cache Ensembl data (default is ~/.cache/denovonear)")
    
    subparsers = parser.add_subparsers()
    
    ############################################################################
    # CLI options for clustering
    cluster = subparsers.add_parser('cluster', parents=[parent],
        description="Tests the proximity of de novo mutations in genes.")
    cluster.add_argument("--in", dest="input", required=True, help="Path to "
        "file listing known mutations in genes. See example file in data folder "
        "for format.")
    
    cluster.set_defaults(func=clustering)
    
    ############################################################################
    # CLI options for identifing transcripts to use
    transcripts = subparsers.add_parser('transcripts', parents=[parent],
        description="Identify transcripts for a gene containing de novo events.")
    transcripts.add_argument("--de-novos", required=True, help="Path to "
        "file listing de novo variants in genes.")
    transcripts.set_defaults(func=find_transcripts)
    group = transcripts.add_mutually_exclusive_group(required=True)
    group.add_argument("--all-transcripts", action="store_true", default=False,
        help="Flag if you want to identify all transcripts with more than "
        "one de novo on it.")
    group.add_argument("--minimise-transcripts", action="store_true",
        default=False, help="Flag if you want to identify the minimal set of "
        "transcripts to contain all de novos.")
    
    ############################################################################
    # CLI options for identifing transcripts to use
    rater = subparsers.add_parser("rates", parents=[parent],
        description="determine mutation rates for genes given transcript IDs.")
    rater.add_argument("--genes", help="Path to file "
        "listing HGNC symbols, with one or more transcript IDs per gene. "
        "The tab-separated input format is gene symbol followed by transcript "
        "ID. Alternative transcripts are listed on separate lines.")
    rater.set_defaults(func=gene_rates)
    
    args = parser.parse_args()
    if 'func' not in args:
        print('Use one of the subcommands: cluster, rates, or transcripts\n')
        parser.print_help()
        sys.exit()
    
    return args

def main():
    
    args = get_options()
    
    ensembl = EnsemblRequest(args.cache_folder, args.genome_build.lower())
    mut_dict = load_mutation_rates(args.rates)
    output = open(args.out, "wt")
    
    args.func(ensembl, mut_dict, output, args)

if __name__ == '__main__':
    main()
