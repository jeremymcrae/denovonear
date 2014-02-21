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
from get_transcript_sequences import GetTranscriptSequence
from weighted_choice import WeightedChoice
from load_mutation_rates import load_trincleotide_mutation_rates
from load_known_de_novos import load_known_de_novos

# define some paths and config files
USER_FOLDER = "/nfs/users/nfs_j/jm33/"
APP_FOLDER = os.path.join(USER_FOLDER, "apps", "mutation_rates")
DATA_FOLDER = os.path.join(APP_FOLDER, "data")

ENSEMBL_TO_HGNC_FILE = os.path.join(DATA_FOLDER, "ensembl_transcript_id_to_hgnc_id.txt")
MUTATION_RATES_FILE = os.path.join(DATA_FOLDER, "forSanger_1KG_mutation_rate_table.txt")
KNOWN_MUTATIONS_FILE = os.path.join(DATA_FOLDER, "DNG_Variants_28Jan2014.xlsx")


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

def identify_transcript(ensembl, transcript_ids):
    """ for a given HGNC ID, finds the transcript with the longest CDS
    """
    
    max_length = 0
    max_transcript_id = None
    
    for transcript_id in transcript_ids:
        # get the transcript's protein sequence via the ensembl REST API
        seq = ensembl.get_protein_seq_for_transcript(transcript_id)
        
        # ignore transcripts without protein sequence
        if seq == "Sequence unavailable":
            continue
        
        # only swap to using the transcript if it is the longest
        if len(seq) > max_length:
            max_length = len(seq)
            max_transcript_id = transcript_id
    
    return max_transcript_id

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
    
    return float(sum(distances))/len(distances)

def build_weighted_site_rates_for_gene(gene, mut_dict):
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
        
        # get the distances to the closest exon boundaries
        exon_start, exon_end = gene.find_closest_exon(bp)
        exon_start_dist = abs(exon_start - bp)
        exon_end_dist = abs(exon_end - bp)
        
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
            # only include the site specific probability if it mutates an
            # amino acid, or occurs close to an intron/exon boundary.
            if current_aa != initial_aa or exon_start_dist < 3 or exon_end_dist < 3:
                mutation_possibilities.append([cds_pos, mut_dict[seq][mutated_seq]])
    
    # throw the mutation probabilities into a random sampler that samples by 
    # the weights of the mutation probabilities
    sampler = WeightedChoice(mutation_possibilities)
    
    return sampler

def build_distance_distribution(gene_weights, sample_n=2, max_iter=100):
    """ creates a distribution of distances between mutations in a single gene
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

def generate_distances_for_gene(hgnc_mapper, ensembl, mut_dict, gene_id, sample_n, max_iter):
    """ sort out all the necessary sequences and mutations rates for a gene
    """
    
    potential_transcript_ids = hgnc_mapper[gene_id]
    transcript_id = identify_transcript(ensembl, potential_transcript_ids)
    
    if transcript_id == None:
        raise ValueError(gene_id + " lacks coding transcripts")
    
    # get the sequence for the identified transcript
    (chrom, start, end, strand, genomic_sequence) = ensembl.get_genomic_seq_for_transcript(transcript_id, expand=10)
    cds_sequence = ensembl.get_cds_seq_for_transcript(transcript_id)
    
    # get the locations of the exons and cds from ensembl
    cds_ranges = ensembl.get_cds_ranges_for_transcript(transcript_id)
    exon_ranges = ensembl.get_exon_ranges_for_transcript(transcript_id)
    
    # start an interval object with the locations and sequence
    transcript = Interval(transcript_id, start, end, strand, chrom, exon_ranges, cds_ranges)
    transcript.add_cds_sequence(cds_sequence)
    transcript.add_genomic_sequence(genomic_sequence, offset=10)
    
    # build the site weights, and return a list of mean dist
    site_weights = build_weighted_site_rates_for_gene(transcript, mut_dict)
    random_distances = build_distance_distribution(site_weights, sample_n, max_iter)
    
    return random_distances, transcript

def main():
    # load all the data
    hgnc_mapper = index_hgnc_to_ensembl(ENSEMBL_TO_HGNC_FILE)
    ensembl = GetTranscriptSequence()
    mut_dict = load_trincleotide_mutation_rates(MUTATION_RATES_FILE)
    
    known_de_novos = load_known_de_novos(KNOWN_MUTATIONS_FILE)
    
    # gene_id = "OR5A1"
    # gene_id = "CTC1"
    # de_novo_events = [1,50]
    # gene_distances = generate_distances_for_gene(hgnc_mapper, ensembl, mut_dict, gene_id, len(de_novo_events), 1000000)
    # gene_distances = generate_distances_for_gene(hgnc_mapper, ensembl, mut_dict, gene_id, len(de_novo_events), 1000000)
    
    tested = 0
    for gene_id in known_de_novos:
        
        de_novo_events = known_de_novos[gene_id]
        # don't analyse genes with only one de novo functional mutation
        if len(de_novo_events) < 2:
            continue
        tested += 1
        
        random_distances, interval = generate_distances_for_gene(hgnc_mapper, ensembl, mut_dict, gene_id, len(de_novo_events), 1000000)
        dist = sorted(random_distances)
        
        # need to convert the de novo event positions into CDS positions
        cds_positions = []
        for pos in de_novo_events:
            try:
                cds_positions.append(interval.get_coding_distance(interval.get_cds_start(), pos))
            except AssertionError:
                (start, end) = interval.find_closest_exon(pos)
                
                start_dist = abs(start - pos)
                end_dist = abs(end - pos)
                
                # if the var is outside the exon, but affects a splice site, 
                # swap it to using the splice site location
                if start_dist < 3:
                    cds_positions.append(interval.get_coding_distance(interval.get_cds_start(), start))
                elif end_dist < 3:
                    cds_positions.append(interval.get_coding_distance(interval.get_cds_start(), end))
                else:
                    raise ValueError
        
        observed_distance = get_mean_distance_between_positions(cds_positions)
        
        pos = bisect.bisect_right(dist, observed_distance)
        print("{0} (n={1})\tdenovo dist: {2:0.1f}\tsim rank: {3}\tn sims: {4}\tsims with dist < observed: {5}".format(gene_id, len(de_novo_events), observed_distance, pos, len(dist), pos/len(dist)))
        
        # # start off by only testing a few genes, so we can be sure it's working 
        # # well
        # if tested > 8:
        #     break
        

if __name__ == '__main__':
    main()





