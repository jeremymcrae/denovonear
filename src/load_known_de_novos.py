""" loads a file containing known de novo mutations
"""

from __future__ import print_function
from __future__ import division

from xlrd import open_workbook

# functional consequences, minus any indel consequences
functional_consequences = ["splice_donor_variant", "splice_acceptor_variant",\
    "initiator_codon_variant", "missense_variant", "transcript_amplification", \
    "stop_gained", "stop_lost", "coding_sequence_variant"]

missense_consequences = ["initiator_codon_variant", "missense_variant",\
    "stop_lost"]

nonsense_consequences = ["splice_donor_variant", "splice_acceptor_variant",\
    "stop_gained"]

def load_known_de_novos(filename):
    """ load known mutations into dict indexed by HGNC ID.
    """
    
    # sometimes the input is an excel file
    if filename.endswith("xlsx"):
        # open an xlsx file containing the current list of de novo mutations in
        # probands studied to date
        book = open_workbook(filename, on_demand=True)
        sheet = book.sheet_by_name("DNMs_nonred")
        header = table.row_values(0)
        
        # convert the excel file to a list of lists
        table = []
        row = 0
        while row < table.nrows - 1:
            row += 1
            values = table.row_values(row)
            table.append(values)
    # otherwise assume the file is a tab-separated text file
    else:
        f = open(filename, "r")
        read = f.readlines()
        header = read[0].strip().split("\t")
        
        # convert the tab separated file to a list of lists
        table = []
        for line in read[1:]:
            table.append(line.rstrip().split("\t"))
    
    # get the positions of the columns that we are interested in
    gene_name_col = header.index("gene_name")
    position_col = header.index("pos")
    consequence_col = header.index("consequence")
    snp_or_indel_col = header.index("snp_or_indel")
    unused_consequences = set([])
    
    genes_dict = {}
    for values in table:
        # print(values)
        gene = values[gene_name_col]
        position = values[position_col]
        consequence = values[consequence_col]
        snp_or_indel = values[snp_or_indel_col]
        
        # don't include de novos that aren't functionally important
        if consequence not in functional_consequences:
            unused_consequences.add(consequence)
            continue
        
        # ignore indels (some splice_acceptor_variants (in the
        # functional_consequences) are indels
        if "INDEL" in snp_or_indel:
            continue
        
        # trim out variants that are missing data
        if gene == "" or gene == "." or position == "NA":
            continue
        
        if gene not in genes_dict:
            genes_dict[gene] = {"functional": [], "missense": [], "nonsense": []}
        
        genes_dict[gene]["functional"].append(int(position))
        if consequence in missense_consequences:
            genes_dict[gene]["missense"].append(int(position))
        elif consequence in nonsense_consequences:
            genes_dict[gene]["nonsense"].append(int(position))
    
    return genes_dict
        
        
        
    
    
    
    