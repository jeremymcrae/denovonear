""" loads a file containing known de novo mutations
"""

from __future__ import print_function
from __future__ import division

# functional consequences, minus any indel consequences
functional_consequences = ["splice_donor_variant", "splice_acceptor_variant",\
    "initiator_codon_variant", "missense_variant", "transcript_amplification", \
    "stop_gained", "stop_lost", "coding_sequence_variant"]

missense_consequences = ["initiator_codon_variant", "missense_variant",\
    "stop_lost"]

nonsense_consequences = ["splice_donor_variant", "splice_acceptor_variant",\
    "stop_gained"]
    
synonymous_consequences = ["synonymous_variant"]

def load_known_de_novos(filename):
    """ load known mutations into dict indexed by HGNC ID.
    """
    
    f = open(filename, "r")
    header = f.readline().strip().split("\t")
    
    # get the positions of the columns that we are interested in
    gene_name_col = header.index("gene_name")
    position_col = header.index("pos")
    consequence_col = header.index("consequence")
    snp_or_indel_col = header.index("snp_or_indel")
    
    genes_dict = {}
    for line in f:
        line = line.rstrip().split("\t")
        gene = line[gene_name_col]
        position = line[position_col]
        consequence = line[consequence_col]
        snp_or_indel = line[snp_or_indel_col]
        
        # ignore indels (some splice_acceptor_variants (in the
        # functional_consequences) are indels
        if "INDEL" in snp_or_indel:
            continue
        
        # trim out variants that are missing data
        if gene == "" or gene == "." or position == "NA":
            continue
        
        if gene not in genes_dict:
            genes_dict[gene] = {"functional": [], "missense": [], \
                "nonsense": [], "synonymous": []}
        
        genes_dict[gene]["functional"].append(int(position))
        if consequence in missense_consequences:
            genes_dict[gene]["missense"].append(int(position))
        elif consequence in nonsense_consequences:
            genes_dict[gene]["nonsense"].append(int(position))
        elif consequence in synonymous_consequences:
            genes_dict[gene]["synonymous"].append(int(position))
        
    return genes_dict
        
        
        
    
    
    
    