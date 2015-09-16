""" script to obtain 1000 genomes variation data in CDS exons
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import gzip
import copy

try:
    import xml.etree.cElementTree as et
except ImportError:
    import xml.etree.ElementTree as et

LOF = {"STOP-GAIN": "stop_gained",
       "Splice Site": "splice_donor_variant",
       "exon_loss": "transcript_ablation",
       "splice-3": "splice_acceptor_variant",
       "splice-5": "splice_donor_variant",
       "splice_acceptor_variant": "splice_acceptor_variant",
       "splice_donor_variant": "splice_donor_variant",
       "stop_gained": "stop_gained"}

missense = {"missense": "missense_variant",
            "STOP-LOSS": "stop_lost",
            "missense_variant": "missense_variant",
            "stop_lost": "stop_lost"}

functional = copy.deepcopy(LOF)
functional.update(missense)

ddg2p_file = "/nfs/ddd0/Data/datafreeze/1133trios_20131218/DDG2P_with_genomic_coordinates_20131107_updated_TTN.tsv"

def open_ddg2p(filename):
    """ get a set of DDG2P genes suitable for benchmarking the clustering code
    
    Args:
        filename: path to DDG2P tab-separated file
    
    Returns:
        set of gene symbols
    """
    
    ddg2p = open(filename, "r")
    
    allowed_mechs = set(["Loss of function", "Dominant negative", "Activating"])
    
    genes = set([])
    for line in ddg2p:
        line = line.split("\t")
        
        gene = line[3]
        gene_type = line[4]
        mode = line[5]
        mech = line[6]
        
        # ignore all the low confidence genes
        if gene_type == "Possible DD Gene":
            continue
        
        # we only want Dominant genes
        if mode != "Monoallelic":
            continue
        
        # and we only want certain mechanisms
        if mech not in allowed_mechs:
            continue
        
        genes.add(gene)
    
    return genes

class ParseClinVar():
    """ parses ClinVar XML data source
    """
    
    def __init__(self, ddg2p, filename="/nfs/users/nfs_j/jm33/apps/denovonear/data/benchmarking/ClinVarFullRelease_2014-04.xml.gz"):
        """ opens a file object for the clinvar file
        """
        
        self.ddg2p_genes = ddg2p
        self.xml_file = gzip.open(filename)
        self.vars = []
    
    def parse_xml(self):
        """ parses the xml file, only looking into ClinVarSet nodes
        """
        
        self.genes = set([])
        
        for (event, elem) in et.iterparse(self.xml_file):
            if elem.tag == "ClinVarSet":
                # print(elem.tag)
                self.parse_clinvarset(elem)
                elem.clear()
        
    def parse_clinvarset(self, clin_elem):
        """ parse a ClinVar set node
        """
        
        clin_var = clin_elem.find("ReferenceClinVarAssertion")
        measure_set = clin_var.find("MeasureSet")
        measure = measure_set.find("Measure")
        
        # ignore non SNVs
        var_type = measure.get("Type")
        if var_type != "single nucleotide variant":
            return
        
        consequence = None
        for attribute_set in measure.findall("AttributeSet"):
            if attribute_set[0].get("Type") == "MolecularConsequence":
                consequence = attribute_set[0].text
        
        measure_rel = measure.find("MeasureRelationship")
        if measure_rel is None:
            measure_rel = measure
        clin_var_ID = clin_elem.get("ID")
        
        (chrom, start, stop, assembly) = self.get_location(measure, measure_rel)
        
        hgnc = None
        symbol = measure_rel.find("Symbol")
        if symbol is not None:
            gene_symbol = symbol.find("ElementValue")
            hgnc = gene_symbol.text
            # print(hgnc)
        
        if start is not None and abs(int(start) - int(stop)) == 0 and \
                consequence in functional and hgnc is not None and hgnc in self.ddg2p_genes:
            self.genes.add(hgnc)
            consequence = functional[consequence]
            self.vars.append((hgnc, chrom, start, consequence, clin_var_ID, assembly))
            # print(hgnc, chrom, start, consequence)
        
        # if hgnc == "ALMS1":
        #     print(clin_var_ID, chrom, start, stop, assembly, hgnc, consequence)
    
    def get_location(self, measure, measure_rel):
        """ parses location data in the xml tree.
        """
        
        chrom = None
        start = None
        stop = None
        assembly = None
        if measure.find("SequenceLocation") is not None:
            locations = measure.findall("SequenceLocation")
            for seq_location in locations:
                if seq_location.get("Assembly") == "GRCh37":
                    chrom = seq_location.get("Chr")
                    start = seq_location.get("start")
                    stop = seq_location.get("stop")
                    assembly = seq_location.get("Assembly")
        
        if start is not None and stop is not None and abs(int(start) - int(stop)) < 2:
            return (chrom, start, stop, assembly)
        
        seq_location = measure_rel.find("SequenceLocation")
        
        chrom = None
        start = None
        stop = None
        assembly = None
        if seq_location is not None:
            chrom = seq_location.get("Chr")
            start = seq_location.get("start")
            stop = seq_location.get("stop")
            assembly = seq_location.get("Assembly")
        
        return (chrom, start, stop, assembly)


def main():
    ddg2p = open_ddg2p(ddg2p_file)
    
    parser = ParseClinVar(ddg2p)
    parser.parse_xml()
    # print(len(parser.genes))
    
    output = open("clinvar_variants.txt", "w")
    output.write("gene_name\tchr\tpos\tconsequence\tsnp_or_indel\tclin_var_ID\tassembly\n")
    for var in sorted(parser.vars):
        line = "{0}\t{1}\t{2}\t{3}\tDENOVO-SNP\t{4}\t{5}\n".format(*var)
        output.write(line)
        
    output.close()

if __name__ == '__main__':
    main()
