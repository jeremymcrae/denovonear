""" loads a file containing known de novo mutations
"""

from __future__ import print_function
from __future__ import division

missense = ["missense_variant", "stop_lost", "inframe_deletion",
    "inframe_insertion", "coding_sequence_variant", "protein_altering_variant"]

lof = ["stop_gained", "splice_acceptor_variant",
    "splice_donor_variant", "frameshift_variant", "initiator_codon_variant",
    "start_lost", "conserved_exon_terminus_variant"]
    
synonymous = ["synonymous_variant"]

def load_de_novos(filename, exclude_indels=True, exclude_snvs=False):
    """ load mutations into dict indexed by HGNC ID.
    """
    
    f = open(filename, "r")
    header = f.readline().strip().split("\t")
    
    genes = {}
    for line in f:
        # ignore header lines
        if line.startswith("gene") or line.lower().startswith("HGNC"):
            continue
        
        line = line.rstrip().split("\t")
        gene = line[0]
        position = int(line[2]) - 1
        consequence = line[3]
        var_type = line[4]
        
        # ignore indels (some splice_acceptor_variants (in the
        # functional_consequences) are indels
        if exclude_indels and "INDEL" in var_type.upper():
            continue
        if exclude_snvs and "SNV" in var_type.upper():
            continue
        
        # trim out variants that are missing data
        if gene == "" or gene == "." or position == "NA":
            continue
        
        if gene not in genes:
            genes[gene] = {"missense": [], "nonsense": []}
        
        if consequence in missense:
            genes[gene]["missense"].append(position)
        elif consequence in lof:
            genes[gene]["nonsense"].append(position)
        
    return genes
        
        
        
    
    
    
    
