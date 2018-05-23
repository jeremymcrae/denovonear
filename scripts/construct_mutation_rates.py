""" Script to generate mutation rates based on local sequence context rates
for Ensembl transcript IDs.
"""

import os
import argparse

from denovonear.load_gene import construct_gene_object
from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.site_specific_rates import SiteRates
from denovonear.frameshift_rate import include_frameshift_rates
from denovonear.log_transform_rates import log_transform

def get_options():
    """ get the command line switches
    """
    
    parser = argparse.ArgumentParser(description="determine mutation rates "
        "for genes given transcript IDs.")
    parser.add_argument("--genes", help="Path to file "
        "listing HGNC symbols, with one or more transcript IDs per gene. "
        "The tab-separated input format is gene symbol followed by transcript "
        "ID. Alternative transcripts are listed on separate lines.")
    
    parser.add_argument("--out", required=True, help="output filename")
    parser.add_argument("--rates",
        help="Path to file containing sequence-context based mutation rates.")
    parser.add_argument("--genome-build", choices=["grch37", "GRCh37",
        "grch38", "GRCh38"], default="grch37", help="Genome build "
        "that the de novo coordinates are based on (currently GRCh37 or GRCh38).")
    parser.add_argument("--cache-folder",
        default=os.path.join(os.path.dirname(__file__), "cache"), help="folder "
            "to cache Ensembl data into (defaults to clustering code directory).")
    
    return parser.parse_args()

def load_genes(path):
    """ load a file listing gene and transcript IDs
    
    Eeach gene can have one or more transcript IDs associated with it, so we
    build a dictionary, indexed by HGNC symbols, and for each gene entry, retain
     a list of the possible transcript IDs.
    
    Args:
        path: path to file containing gene and transcript IDs, with each unique
            transcript ID for a gene on different lines eg
            
            gene_1    transcript_1.1    length_1    denovo_count
            gene_1    transcript_1.2    length_2    denovo_count
            gene_1    transcript_1.2    length_3    denovo_count
            gene_2    transcript_2.1    length_3    denovo_count
    
    Returns:
        dict of transcripts eg {'CTC1': ["ENST00000315684", "ENST00000485511"]}
    """
    
    transcripts = {}
    with open(path, "r") as f:
        for line in f:
            if line.startswith("hgnc"):
                continue
            
            symbol, tx_id, _, _ = line.strip().split("\t")
            
            if symbol not in transcripts:
                transcripts[symbol] = []
            
            transcripts[symbol].append(tx_id)
            
    return transcripts

def get_mutation_rates(gene_id, transcripts, mut_dict, ensembl):
    """ determines missense, nonsense and synonymous mutation rates for a gene
    
    This can estimate a mutation rate from the union of transcripts for a gene.
    This is a biased estimate of the mutation rate, where the mutation rate
    estimates is biased towards the rate from the first-ranked transcripts,
    which I prioritise by how many de novos they contain, and how long the
    coding sequence is.
    
    This isn't a problem when different transcripts have the same coding
    sequence within their shared regions, as the rates will come outthe same,
    but may differ two transcript share an overlapping region, but not in the
    same frame, so that the sites that are missense, and nonsense will differ
    between transcripts, and thus would produce different estimates of the
    mutation rate.
    
    Args:
        gene_id: ID for the current gene (can be a transcript ID, if we are
            examining single transcripts only, or can be a HGNC ID, if we are
            examining the union of mutation rates from multiple transcripts for
            a single gene).
        transcripts: dictionary of transcripts for a gene, indexed by gene_id
        mut_dict: dictionary of local sequence context mutation rates
        ensembl: EnsemblRequest object, to retrieve information from Ensembl.
    
    Returns:
        tuple of (missense, nonsense, synonymous) mutation rates
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

def main():
    
    args = get_options()
    
    # load all the data
    ensembl = EnsemblRequest(args.cache_folder, args.genome_build)
    mut_dict = load_mutation_rates(args.rates)
    
    transcripts = load_genes(args.genes)
    
    output = open(args.out, "w")
    header = ['transcript_id', 'chrom', 'length', 'missense_rate', 'nonsense_rate',
        'splice_lof_rate', 'splice_region_rate', 'synonymous_rate']
    output.write('\t'.join(header) + "\n")
    
    for symbol in sorted(transcripts):
        print(symbol)
        try:
            rates, tx, length = get_mutation_rates(symbol, transcripts[symbol],
                mut_dict, ensembl)
            # log transform rates, for consistency with Samocha et al.
            line = "{}\t{}\t{}\t{}".format(symbol, tx.get_chrom(), length, log_transform(rates))
        except (ValueError, KeyError) as error:
            print("{}\t{}\n".format(symbol, error))
            line = "{}\t{}\tNA\tNA\tNA\tNA\tNA\tNA".format(symbol, tx.get_chrom())
        
        output.write(line + '\n')
    
    output.close()
    include_frameshift_rates(args.out)

if __name__ == '__main__':
    main()
