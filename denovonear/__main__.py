""" Script to investigate the probability of multiple mutations clustering
within a single gene.
"""

import os
import sys
import asyncio
import argparse
import logging

from denovonear.rate_limiter import RateLimiter
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.load_de_novos import load_de_novos
from denovonear.cluster_test import cluster_de_novos

from denovonear.load_gene import (load_gene, construct_gene_object,
    count_de_novos_per_transcript, minimise_transcripts)
from denovonear.site_specific_rates import SiteRates
from denovonear.frameshift_rate import include_frameshift_rates
from denovonear.log_transform_rates import log_transform
from denovonear.gencode import Gencode

async def _load_gencode(symbols):
    ''' load gene coords and sequence via ensembl
    '''
    gencode = Gencode()
    async with RateLimiter(per_second=15) as ensembl:
        tasks = [load_gene(ensembl, symbol) for symbol in symbols]
        genes = await asyncio.gather(*tasks)
        for gene in genes:
            gencode.add_gene(gene)
        return gencode

def load_gencode(symbols, gencode=None, fasta=None):
    ''' load genes from gencode annotations file, with ensembl as backup
    '''
    if gencode and fasta:
        return Gencode(gencode, fasta)
    
    # use ensembl as backup if gencode file not available. This restricts the
    # asynchronous calls to within one section, and ensures they are called together
    return asyncio.get_event_loop().run_until_complete(_load_gencode(symbols))

def clustering(args, output):
    
    mut_dict = load_mutation_rates(args.rates)
    de_novos = load_de_novos(args.input)
    gencode = load_gencode(de_novos, args.gencode, args.fasta)
    
    output.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")
    
    iterations = 1000000
    for symbol in sorted(de_novos):
        logging.info(f'checking {symbol}')
        if len(de_novos[symbol]["missense"] + de_novos[symbol]["nonsense"]) < 2:
            continue
        
        if symbol not in gencode:
            logging.info(f'cannot find {symbol} in gencode')
            continue
        
        probs = cluster_de_novos(symbol, de_novos[symbol], gencode[symbol], iterations, mut_dict)
        
        if probs is None:
            continue
        
        output.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense",
            len(de_novos[symbol]["missense"]), probs["miss_dist"], probs["miss_prob"]))
        output.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "nonsense",
            len(de_novos[symbol]["nonsense"]), probs["nons_dist"], probs["nons_prob"]))

def find_transcripts(args, output):
    
    mut_dict = load_mutation_rates(args.rates)
    de_novos = load_de_novos(args.de_novos)
    gencode = load_gencode(de_novos, args.gencode, args.fasta)
    
    output.write("hgnc_symbol\ttranscript_id\tlength\tde_novos\n")
    
    for symbol in sorted(de_novos):
        logging.info(f'checking {symbol}')
        func_events = de_novos[symbol]["missense"] + de_novos[symbol]["nonsense"]
        
        transcripts = gencode[symbol].transcripts
        
        # find the counts per transcript, depending on whether we want to count
        # for all transcripts containing one or more de novos, or to find the
        # minimum set of transcripts to contain the de novos
        try:
            if args.all_transcripts:
                counts = count_de_novos_per_transcript(transcripts, func_events)
            elif args.minimal_transcripts:
                counts = minimise_transcripts(transcripts, func_events)
        except (ValueError, IndexError):
            logging.error(f"error occured with {symbol}")
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

def get_mutation_rates(transcripts, mut_dict):
    """ determines mutation rates per functional category for transcripts
    
    Args:
        transcripts: list of transcript IDs for a gene
        mut_dict: dictionary of local sequence context mutation rates
    
    Returns:
        tuple of (rates, merged transcript, and transcript CDS length)
    """
    
    rates = {'missense': 0, 'nonsense': 0, 'splice_lof': 0,
        'splice_region': 0, 'synonymous': 0}
    combined = None
    
    for tx in transcripts:
        
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
    
    length = combined.get_coding_distance(combined.get_cds_end())['pos']
    
    return rates, combined, length

def gene_rates(args, output):
    ''' calculate per-consequence mutation rates per gene
    '''
    mut_dict = load_mutation_rates(args.rates)
    transcripts = load_genes(args.genes)
    gencode = load_gencode(transcripts, args.gencode, args.fasta)
    
    header = ['transcript_id', 'chrom', 'length', 'missense_rate', 'nonsense_rate',
        'splice_lof_rate', 'splice_region_rate', 'synonymous_rate']
    output.write('\t'.join(header) + "\n")
    
    for symbol in sorted(transcripts):
        tx_ids = set(transcripts[symbol])
        gene = gencode[symbol]
        txs = [x for x in gene.transcripts if x.get_name().split('.')[0] in tx_ids]
        try:
            rates, tx, length = get_mutation_rates(txs, mut_dict)
            # log transform rates, for consistency with Samocha et al.
            line = "{}\t{}\t{}\t{}".format(symbol, tx.get_chrom(), length, log_transform(rates))
        except (ValueError, KeyError) as error:
            logging.error(f"{symbol}\t{error}\n")
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
    parent.add_argument("--out", default=sys.stdout, help="output filename")
    parent.add_argument("--rates",
        help="optional path to file containing sequence context-based mutation rates.")
    parent.add_argument("--gencode",
        help="optional path to gencode annotations file. If not provided, gene " \
            "coordinates will be obtained via the Ensembl REST API.")
    parent.add_argument("--fasta",
        help="optional path to genome fasta file. If not provided, gene " \
            "coordinates will be obtained via the Ensembl REST API.")
    parent.add_argument("--genome-build", choices=["grch37",
        "GRCh37", "grch38", "GRCh38"], default="grch37", help="Genome build "
        "that the de novo coordinates are based on (GRCh37 or GRCh38")
    parent.add_argument("--log", default=sys.stdout, help="where to write log files")
    
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
    # CLI options for getting mutation rates per gene
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

def open_output(path):
    ''' open output, which could be standard out
    '''
    try:
        output = open(path, 'wt')
    except TypeError:
        output = path
    return output

def main():
    args = get_options()
    FORMAT = '%(asctime)-15s %(message)s'
    log = open(args.log, 'at') if isinstance(args.log, str) else args.log
    logging.basicConfig(stream=log, format=FORMAT, level=logging.INFO)
    
    output = open_output(args.out)
    args.func(args, output)

if __name__ == '__main__':
    main()
