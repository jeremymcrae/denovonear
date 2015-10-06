""" loads a file containing known de novo mutations
"""

from __future__ import print_function
from __future__ import division

missense_consequences = ["initiator_codon_variant", "missense_variant",\
    "stop_lost", "inframe_deletion", "inframe_insertion"]

nonsense_consequences = ["splice_donor_variant", "splice_acceptor_variant",\
    "stop_gained", "frameshift_variant"]
    
synonymous_consequences = ["synonymous_variant"]

def load_de_novos(filename, exclude_indels=True, exclude_snvs=False):
    """ load mutations into dict indexed by HGNC ID.
    """
    
    f = open(filename, "r")
    header = f.readline().strip().split("\t")
    
    genes_dict = {}
    for line in f:
        # ignore header lines
        if line.startswith("gene") or line.lower().startswith("HGNC"):
            continue
        
        line = line.rstrip().split("\t")
        gene = line[0]
        chrom = line[1]
        position = int(line[2]) - 1
        consequence = line[3]
        snp_or_indel = line[4]
        
        # ignore indels (some splice_acceptor_variants (in the
        # functional_consequences) are indels
        if exclude_indels and "INDEL" in snp_or_indel.upper():
            continue
        if exclude_snvs and "SNV" in snp_or_indel.upper():
            continue
        
        # trim out variants that are missing data
        if gene == "" or gene == "." or position == "NA":
            continue
        
        if gene not in genes_dict:
            genes_dict[gene] = {"functional": [], "missense": [], \
                "nonsense": [], "synonymous": []}
        
        genes_dict[gene]["functional"].append(position)
        if consequence in missense_consequences:
            genes_dict[gene]["missense"].append(position)
        elif consequence in nonsense_consequences:
            genes_dict[gene]["nonsense"].append(position)
        elif consequence in synonymous_consequences:
            genes_dict[gene]["synonymous"].append(position)
        
    return genes_dict
        
        
        
    
    
    
    
