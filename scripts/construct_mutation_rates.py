""" Script to generate mutation rates based on local sequence context rates
for Ensembl transcript IDs.
"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import os
import copy
import math
import argparse
import tempfile

from denovonear.load_gene import construct_gene_object
from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_mutation_rates import load_trincleotide_mutation_rates
from denovonear.site_specific_rates import SiteRates

def get_options():
    """ get the command line switches
    """
    
    parser = argparse.ArgumentParser(description="determine mutation rates \
        for genes given transcript IDs.")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--transcripts", dest="transcript_input",
        help="Path to file listing Ensembl transcript IDs, one ID per line.")
    group.add_argument("--genes", dest="gene_input", help="Path to file" + \
        " listing HGNC symbols, with one or more transcript IDs per gene. " + \
        "The tab-separated input format is gene symbol followed by transcript " + \
        "ID. Alternative transcripts are listed on separate lines.")
    
    parser.add_argument("--out", dest="output", required=True, help="output \
        filename")
    parser.add_argument("--rates", dest="mut_rates", required=True, \
        help="Path to file containing trinucleotide mutation rates.")
    parser.add_argument("--genome-build", dest="genome_build", choices=["grch37",
        "GRCh37", "grch38", "GRCh38"], default="grch37", help="Genome build "+ \
        "that the de novo coordinates are based on (currently GRCh37 or GRCh38).")
    parser.add_argument("--cache-folder", dest="cache_folder", \
        default=os.path.join(os.path.dirname(__file__), "cache"), help="folder \
        to cache Ensembl data into (defaults to clustering code directory).")
    
    parser.add_argument("--coverage-adjust", default=False, action="store_true", \
        help="whether to adjust site mutation rates for sequencing coverage.")
    parser.add_argument("--coverage-dir", help="location of ExAC coverage files")
    
    args = parser.parse_args()
    
    return args.transcript_input, args.gene_input, args.output, args.mut_rates, \
        args.cache_folder, args.genome_build.lower(), args.coverage_adjust, \
        args.coverage_dir

def load_transcripts(path):
    """ load a file listing transcript IDs per line
    
    Args:
        path: path to file containing transcript IDs, one per line
    
    Returns:
        list of transcript IDs eg ["ENST00000315684", "ENST00000485511"]
    """
    
    transcripts = {}
    with open(path, "r") as f:
        for line in f:
            transcripts[line.strip()] = [line.strip()]
            
    return transcripts
    
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
        list of transcript IDs eg ["ENST00000315684", "ENST00000485511"]
    """
    
    transcripts = {}
    with open(path, "r") as f:
        for line in f:
            if line.startswith("hgnc"):
                continue
            
            line = line.strip().split("\t")
            gene_id = line[0]
            transcript_id = line[1]
            
            if gene_id not in transcripts:
                transcripts[gene_id] = []
            
            transcripts[gene_id].append(transcript_id)
            
    return transcripts

def get_mutation_rates(gene_id, transcripts, mut_dict, ensembl, use_cov, cov_dir):
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
    
    missense = 0
    nonsense = 0
    splice_lof = 0
    splice_region = 0
    synonymous = 0
    combined_transcript = None
    
    for transcript_id in transcripts[gene_id]:
        
        # get the gene coordinates, sequence etc, but if the transcript is
        # unusable (hence raises an error), simply move to the next transcript
        try:
            transcript = construct_gene_object(ensembl, transcript_id)
        except ValueError:
            continue
        
        if len(transcript.cds_sequence) % 3 != 0:
            raise ValueError("anomalous_coding_sequence")
        
        # ignore mitochondrial genes, since mitochondiral mutation rates differ
        # from autosomal and allosomal mutation rates
        if transcript.get_chrom() == "MT":
            continue
        
        if combined_transcript is None:
            site_weights = SiteRates(transcript, mut_dict, use_coverage=use_cov)
            combined_transcript = copy.deepcopy(transcript)
        else:
            site_weights = SiteRates(transcript, mut_dict, masked_sites=combined_transcript, use_coverage=use_cov)
            combined_transcript += transcript
        
        if cov_dir is not None:
            site_weights.set_coverage_dir(cov_dir)
        
        missense_rates = site_weights.get_missense_rates_for_gene()
        nonsense_rates = site_weights.get_nonsense_rates_for_gene()
        splice_lof_rates = site_weights.get_splice_loss_of_function_rates_for_gene()
        splice_region_rates = site_weights.get_splice_region_rates_for_gene()
        synonymous_rates = site_weights.get_synonymous_rates_for_gene()
        
        # if any sites have been sampled in the transcript, then add the
        # cumulative probability from those sites to the approporiate
        # mutation rate. Sometimes we won't have any sites for a transcript, as
        # all the sites will have been captured in previous transcripts.
        if len(missense_rates.choices) > 0:
            missense += missense_rates.cum_probs[-1]
        if len(nonsense_rates.choices) > 0:
            nonsense += nonsense_rates.cum_probs[-1]
        if len(splice_lof_rates.choices) > 0:
            splice_lof += splice_lof_rates.cum_probs[-1]
        if len(splice_region_rates.choices) > 0:
            splice_region += splice_region_rates.cum_probs[-1]
        if len(synonymous_rates.choices) > 0:
            synonymous += synonymous_rates.cum_probs[-1]
    
    chrom = combined_transcript.get_chrom()
    length = "NA"
    if combined_transcript is not None:
        length = combined_transcript.get_coding_distance(combined_transcript.cds_min, combined_transcript.cds_max)
    
    return (chrom, length, missense, nonsense, splice_lof, splice_region, synonymous)

def log_transform(values):
    """ log transform a numeric value, unless it is zero, or negative
    """
    
    transformed = []
    for value in values:
        try:
            value = math.log10(value)
        except ValueError:
            value = "NA"
        transformed.append(value)
    
    return transformed

def include_indel_rates(path):
    """ add per-gene indel mutation rates to the output file
    
    We estimate the indel mutation rate as per Nature Genetics 46:944-950
    (2014) doi:10.1038/ng.3050, in which the overall (ie summed across all
    genes) indel mutation rate is estimated as 1.25 x the overall nonsense
    rate. Given this, the indel mutation rate is portioned out to each gene,
    according to the proportion of each gene's CDS length to the length of all
    CDS regions from all genes.
    
    Note that this is an approximation, and doesn't allow for divergence due to
    the base compostion of the CDS.
    
    Args:
        path: path to the output mutation rates (for the nonsense, missense etc)
    """
    
    # we copy the rates file to a temporary file, so that we can stream through
    # it later while writing amended lines to the correct rates path
    temp = tempfile.TemporaryFile(mode="w+b")
    
    # run through the file of rates to find the overall nonsense mutation rate,
    # and the total length of CDS regions in the file.
    nonsense_sum = 0
    length_sum = 0
    with open(path) as handle:
        nonsense_pos = None
        length_pos = None
        for line in handle:
            temp.write(line)
            if line.startswith("transcript_id"):
                nonsense_pos = line.split("\t").index("nonsense_rate")
                length_pos = line.split("\t").index("length")
                continue
            
            line = line.strip().split("\t")
            nonsense_sum += 10**float(line[nonsense_pos])
            length_sum += int(line[length_pos])
    
    handle.close()
    temp.seek(0)
    
    # add the frameshift rates to each line in turn, while writing the output
    # back to the output path
    frameshift_sum = nonsense_sum * 1.25
    with open(path, "w") as handle:
        for line in temp:
            line = line.strip().split("\t")
            if line[0] == "transcript_id":
                line.append("frameshift_rate")
            else:
                # estimate the frameshift rate for the gene
                frameshift_rate = (float(line[2])/length_sum) * frameshift_sum
                frameshift_rate = math.log10(frameshift_rate)
                line.append(str(frameshift_rate))
            
            line = "\t".join(line) +"\n"
            handle.write(line)
    
    temp.close()

def main():
    
    input_transcripts, input_genes, output_file, rates_file, cache_dir, \
        genome_build, use_coverage, coverage_dir = get_options()
    
    # load all the data
    ensembl = EnsemblRequest(cache_dir, genome_build)
    mut_dict = load_trincleotide_mutation_rates(rates_file)
    
    if input_transcripts is not None:
        transcripts = load_transcripts(input_transcripts)
    else:
        transcripts = load_genes(input_genes)
    
    output = open(output_file, "w")
    output.write("transcript_id\tchrom\tlength\tmissense_rate\tnonsense_rate\tsplice_lof_rate\tsplice_region_rate\tsynonymous_rate\n")
    
    for gene_id in sorted(transcripts):
        print(gene_id)
        try:
            rates = get_mutation_rates(gene_id, transcripts, mut_dict, ensembl, \
                use_coverage, coverage_dir)
            
            chrom = rates[0]
            length = rates[1]
            rates = rates[2:]
            # log transform the rates, to keep them consistent with the rates from
            # Daly et al.
            line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(gene_id, chrom, length, *log_transform(rates))
        except ValueError as error:
            line = "{0}\t{1}\n".format(gene_id, error)
        except KeyError as error:
            # ignore genes with odd genomic sequence eg ENST00000436041 in GRCh37
            continue
        
        output.write(line)
    
    output.close()
    
    include_indel_rates(output_file)

if __name__ == '__main__':
    main()
