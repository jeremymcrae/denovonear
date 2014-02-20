""" Script to investigate the probability of multiple mutations clustering 
within a single gene.

simulate multiple mutations within a single gene using mutation rates
    only consider nonsynonymous mutations
    check functional type

simulate clustering of mutations
    sum the probabilities at each variant in the CDS for nonsynonymous mutations
    sum all the probabilities across all the positions in the CDS
    the likelihood of a specific mutation at a particular location is the 
    ratio of that specific trinucleotide change to the summed gene probability

cluster by distance:
    determine distance between two mutations (distance defined as within 
        coding sequence for longest CDS (if > one transcript))

    Rationale is:
        de novos closer together are increasingly unlikely
        multiple de novos hitting same position is near impossible
        three or more de novos - use mean distance between variants

compare clustering of simulated mutations to known mutations
    eg only 1% of the time we obtain de novos within 10 bp of each other

"""


from __future__ import print_function
from __future__ import division
import sys
import os
import bisect
import itertools
import time

from interval import Interval
from prepare_sequence_data import PrepareSequenceData
from weighted_choice import WeightedChoice
from load_mutation_rates import load_trincleotide_mutation_rates

# define some paths and config files
USER_FOLDER = "/nfs/users/nfs_j/jm33/"
APP_FOLDER = os.path.join(USER_FOLDER, "apps", "mutation_rates")
DATA_FOLDER = os.path.join(APP_FOLDER, "data")

ENSEMBL_TO_HGNC_FILE = os.path.join(DATA_FOLDER, "ensembl_transcript_id_to_hgnc_id.txt")
# TRANSCRIPTS_FILE = os.path.join(DATA_FOLDER, "gencode.v19.annotation.bed")
# PROTEIN_SEQUENCES_FILE = os.path.join(DATA_FOLDER, "ensembl_protein_sequence_by_ensembl_id.fa")
# CDS_SEQUENCES_FILE = os.path.join(DATA_FOLDER, "ensembl_fasta.cds.fa")
# GENOMIC_SEQUENCES_FILE = os.path.join(DATA_FOLDER, "ensembl_fasta.genomic.fa")
MUTATION_RATES_FILE = os.path.join(DATA_FOLDER, "forSanger_1KG_mutation_rate_table.txt")


# def load_transcript_from_bed_file(filename, transcript_id):
#     """ opens a bed file as a dict of Interval objects, indexed by name
#     """
    
#     intervals = {}
#     with open(filename) as f:
#         for line in f:
#             line = line.strip()
#             if transcript_id in line:
#                 return Interval(line)
    
#     return None

def index_hgnc_to_ensembl(filename):
    """ build dict to convert HGNC IDs to all their ensembl transcript IDs
    """
    
    mapper = {}
    with open(filename) as f:
        for line in f:
            line = line.strip().split("\t")
            # ignore lines without HGNC IDs
            if len(line) == 1:
                continue
            
            transcript_id = line[0]
            hgnc_id = line[1]
            
            if hgnc_id not in mapper:
                mapper[hgnc_id] = set([])
                
            mapper[hgnc_id].add(transcript_id)
        
    return mapper

def get_mean_distance_between_positions(positions):
    """ gets the mean distance between two or more CDS positions
    """
    
    if len(positions) == 1:
        raise ValueError("only one position given, cannot calculate distance")
    elif len(positions) == 2:
        return abs(positions[0] - positions[1])
    
    pos_pairs = itertools.combinations(positions, 2)
    
    distances = []
    for pos_1, pos_2 in pos_pairs:
        distance = abs(pos_1 - pos_2)
        distances.append(distance)
    
    return float(sum(distances)/len(distances))

def build_weighted_site_rates_for_gene(gene):
    """ build a list of sites in a gene that can give missense mutations, along 
    with their weighted probability of the mutation occuring.
    """
    
    mutation_possibilities = []
    for bp in range(gene.get_cds_start(), gene.get_cds_end()):
        if not gene.in_coding_region(bp):
            continue
            
        bases = {"A", "C", "G", "T"}
        
        cds_pos = gene.convert_genomic_position_to_relative_to_ATG_in_cds(bp)
        if gene.strand == "+":
            cds_pos += 1
        
        codon_number = gene.get_codon_number_for_cds_position(cds_pos)
        pos_within_codon = gene.get_position_within_codon(cds_pos)
        codon = gene.get_codon_sequence(codon_number)
        
        seq = gene.get_trinucleotide_around_cds_position(cds_pos)
        initial_nucleotide = seq[1]
        
        # drop the initial nucleotide, since we want to mutate to other bases
        bases.remove(initial_nucleotide)
        
        initial_aa = gene.translate_codon(codon)
        for base in sorted(bases):
            mutated_seq = seq[0] + base + seq[2]
            
            new_codon = list(codon)
            new_codon[pos_within_codon] = base
            new_codon = "".join(new_codon)
            
            current_aa = gene.translate_codon(new_codon)
            if current_aa != initial_aa:
                mutation_possibilities.append([cds_pos, mut_dict[seq][mutated_seq]])
    
    # throw the mutation probabilities into a random sampler that samples by 
    # the weights of the mutation probabilities
    sampler = WeightedChoice(mutation_possibilities)
    
    return sampler

def build_distance_distribution(gene_weights, sample_n=2, max_iter=100):
    """ creates a distribution of 
    """
    
    distribution = []
    iteration = 0
    while iteration < max_iter:
        iteration += 1
        
        positions = []
        while len(positions) < sample_n:
            site = gene_weights.choice()
            positions.append(site)
        
        distance = get_mean_distance_between_positions(positions)
        distribution.append(distance)
    
    return distribution


def generate_distances_for_gene(sequences, mut_dict, gene_id, sample_n, max_iter):
    """ sort out all the necessary sequences and mutations rates for a gene
    """
    
    transcript_id = sequences.identify_transcript(gene_id)
    
    if transcript_id == None:
        raise ValueError(gene_id + " lacks coding transcripts")
    
    # transcript = load_transcript_from_bed_file(transcript_file, transcript_id)
    genomic_sequence = sequences.get_genomic_sequence(transcript_id, expand=10)
    cds_sequence = sequences.get_cds_sequence(transcript_id)
    
    (start, end, strand, chrom) = sequences.get_gene_start_and_end_for_transcript(transcript_id, gene_id)
    cds_ranges = sequences.get_cds_ranges_for_transcript(transcript_id)
    exon_ranges = sequences.get_exon_ranges_for_transcript(transcript_id)
    
    transcript = Interval(transcript_id, start, end, strand, chrom, exon_ranges, cds_ranges)
    
    transcript.add_cds_sequence(cds_sequence)
    transcript.add_genomic_sequence(genomic_sequence, offset=10)
    
    # build the site weights, and return a list of mean dist
    site_weights = build_weighted_site_rates_for_gene(transcript)
    random_distances = build_distance_distribution(site_weights, sample_n, max_iter)
    
    return random_distances


# load all the data
hgnc_mapper = index_hgnc_to_ensembl(ENSEMBL_TO_HGNC_FILE)
sequences = PrepareSequenceData(hgnc_mapper)
mut_dict = load_trincleotide_mutation_rates(MUTATION_RATES_FILE)


gene_id = "TCOF1"
observed_distance = 50
# gene_id = "OR5A1"

gene_distances = generate_distances_for_gene(sequences, mut_dict, gene_id, 3, 1000000)

dist = sorted(gene_distances)

pos = bisect.bisect_left(dist, observed_distance)
print(pos, len(dist), pos/len(dist))




