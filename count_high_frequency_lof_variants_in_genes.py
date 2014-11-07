""" counts high-frequency LOF variants in genes in human genome.

Uses 1000 Genomes variation data to find variants in genes, and estimate the 
minor allele frequency in continental populations.
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
import os
import json
import logging
import argparse

from src.ensembl_requester import EnsemblRequest
from src.load_gene import get_deprecated_gene_ids, get_transcript_lengths, \
    construct_gene_object, get_de_novos_in_transcript, load_gene
from src.ensembl_requester import EnsemblRequest
from src.thousand_genomes.thousand_genomes_variants import Extract1000Genomes

logging.basicConfig(filename='high_freq_lof_variants.log',level=logging.WARNING)

CONSEQUENCES = {"lof": set(["transcript_ablation","splice_donor_variant", \
    "splice_acceptor_variant", "frameshift_variant", "stop_gained", \
    "coding_sequence_variant"]), \
    "missense": set(["initiator_codon_variant", "inframe_insertion", 
    "inframe_deletion", "missense_variant", "transcript_amplification", 
    "stop_lost", "splice_region_variant"]), \
    "synonymous": set(["synonymous_variant"])}

def get_options():
    """
    """
    
    parser = argparse.ArgumentParser(description="Counts common LoF variants in genes.")
    parser.add_argument("--out", dest="output", required=True, help="output \
        filename")
    parser.add_argument("--min-freq", dest="min_freq", default=0.05, help="Minimum accepted minor genotype frequency")
    parser.add_argument("--max-freq", dest="max_freq", default=0.6, help="Maximum accepted minor genotype frequency")
    parser.add_argument("--hgnc", dest="hgnc_filename", default="/nfs/users/nfs_j/jm33/gel_annotation_exercise/data/ensembl.gene_symbols_and_chrom_positions.txt", help="file listing genes to analyse")
    parser.add_argument("--cache-folder", dest="cache_folder", \
        default=os.path.join(os.path.dirname(__file__), "cache"), help="folder \
        to cache Ensembl data into (defaults to clustering code directory)")
    parser.add_argument("--genome-build", dest="genome_build", choices=["grch37",
        "GRCh37", "grch38", "GRCh38"], default="grch37", help="Genome build "+ \
        "that the de novo coordinates are based on (GrCh37 or GRCh38")
    
    args = parser.parse_args()
    
    return args.output, float(args.min_freq), float(args.max_freq), args.hgnc_filename, args.cache_folder, args.genome_build

class EnsemblWithVariants(EnsemblRequest):
    
    def parse_consequence(self, transcript_id, request_data):
        """ parse json request data from ensembl for consequence in a transcript
        """
        
        parsed_json = json.loads(request_data)[0]
        consequence = parsed_json["most_severe_consequence"]
        
        return consequence
    
    def get_variant_by_pos(self, chrom, start_pos, ref_allele, alt_allele, gene):
        """ gets the consequence for a variant from a sequence location
        """
        
        transcript_id = gene.name
        if gene.strand == "-":
            start_pos += 1
        
        # if the variant is a indel, figure out the end position, so we can get
        # the correct variant consequence
        end_pos = start_pos
        if len(ref_allele) > 1:
            end_pos += len(ref_allele) - 1
        
        # some rare structural CNVs surround the entire gene, but their size
        # is challenging for the REST API to return a value. Instead of trying
        # to request an overlay large variant's consequence, just return a
        # default value.
        if start_pos < gene.start and end_pos > gene.end:
            return "transcript_ablated"
        
        ref_seq = self.get_genomic_seq_for_region(chrom, start_pos, end_pos)
        if ref_seq != ref_allele and ref_seq == alt_allele:
            alt_allele = ref_allele
            ref_allele = ref_seq
        
        headers = {"content-type": "application/json"}
        
        self.request_attempts = 0
        ext = "/vep/human/region/{0}:{1}-{2}:1/{3}".format(chrom, start_pos, end_pos, alt_allele)
        r = self.ensembl_request(ext, "{0}:{1}-{2}".format(chrom, start_pos, end_pos), headers)
        
        consequence = self.parse_consequence(transcript_id, r)
        
        return consequence
    
    def get_variant_by_id(self, var_id, transcript_id, chrom, pos, ref_allele, alt_allele, gene):
        """ gets the consequence for a variant based on a rsID
        """
        
        headers = {"content-type": "application/json"}
        
        self.request_attempts = 0
        ext = "/vep/human/id/{0}/consequences?".format(var_id)
        try:
            r = self.ensembl_request(ext, var_id, headers)
        except ValueError:
            return self.get_variant_by_pos(chrom, pos, ref_allele, alt_allele, gene)
        
        consequence = self.parse_consequence(transcript_id, r)
        
        return consequence

def load_hgnc_symbols(filename):
    """ load a file containing all the current HGNC symbols
    """
    
    genes = set([])
    with open(filename, "r") as f:
        for line in f:
            line = line.split("\t")
            hgnc = line[0]
            
            # ignore blank lines
            if hgnc == "":
                continue
            
            genes.add(hgnc)
    
    return sorted(genes)

def get_already_loaded_genes(filename):
    """find the genes which have already been assessed, so we can resume from 
    where we left off (if the code or system breaks)
    """
    
    prior_genes = set([])
    if os.path.exists(filename):
        current_genes = open(filename, "r")
        for line in current_genes:
            line = line.strip().split("\t")
            gene = line[0]
            prior_genes.add(gene)
        current_genes.close()
    
    return prior_genes

def main():
    
    output_filename, min_freq, max_freq, hgnc_filename, cache_dir, genome_build = get_options()
    
    thousand_genomes = Extract1000Genomes("/lustre/scratch113/teams/hurles/users/jm33/1000genomes/", frequency="genotype")
    ensembl = EnsemblWithVariants(cache_dir, genome_build)
    hgnc_symbols = load_hgnc_symbols(hgnc_filename)
    
    # find the genes which have already been assessed, so we can resume from 
    # where we left off (if the code or system breaks)
    prior_genes = get_already_loaded_genes(output_filename)
    
    output = open(output_filename, "a")
    
    for hgnc in hgnc_symbols:
        # hgnc = "GPR107"
        logging.warning("loading: {0}".format(hgnc))
        if hgnc in prior_genes:
            continue
        prior_genes.add(hgnc)
        
        # load the transcript ID, exon positions and sequence. If we raise an
        # error, just move to the next gene (eg ABO, which tries to load a
        # transcript without protein coding sequence, rather than a transcript 
        # on a patch contig).
        try:
            gene = load_gene(ensembl, hgnc)
        except (IndexError, ValueError) as e:
            logging.warning("unable to load: {0}".format(hgnc))
            continue
        
        thousand_genomes.set_gene(gene)
        func = thousand_genomes.filter_variants(min_maf=min_freq, max_maf=max_freq, ignore_indels=False)
        chrom = gene.chrom
        transcript_id = gene.name
        # for category in func:
        #     for (var, ref, alt) in func[category]:
        #         if var.id != ".":
        #             consequence = ensembl.get_variant_by_id(var.id, transcript_id, chrom, var.pos, ref, alt, gene)
                
        #         if consequence not in CONSEQUENCES[category]:
        #             exon_start, exon_end = gene.find_closest_exon(var.pos)
        #             boundary_dist = min(abs(exon_start - var.pos), abs(exon_end - var.pos))
        #             print(var.id, var.pos, category, consequence, boundary_dist)
        
        # break
        output.write(hgnc + "\t" + str(len(func["lof"])) + "\n")

if __name__ == "__main__":
    main()





