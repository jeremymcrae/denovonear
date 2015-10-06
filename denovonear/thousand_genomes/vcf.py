""" class to interact with thousand genomes VCF data
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os

import pysam

class VCF(object):
    """ intended to be a subclass of Extract1000Genomes, where this contains
    functions to connect to the 1000 Genomes VCF, and extract variant records
    """
    
    def set_gene(self, gene):
        """ set a gene for analysis - pull out CDS variants from the VCF
        
        Args:
            gene: Transcript object, containing gene coordinates and sequence
        """
        
        self.gene = gene
        chrom = self.gene.get_chrom()
        
        if os.path.exists(self.vcf_folder):
            chrom_file = "ALL.chr" + str(chrom) + ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
        else:
            chrom_file = "ALL.chr" + str(chrom) + ".integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
        tabix_filename = self.vcf_folder + chrom_file
        
        # set up a tabix connection to the vcf file, so we can efficiently
        # extract variants in a given region
        vcf_tabix = pysam.VCF()
        vcf_tabix.connect(tabix_filename)
        
        # extract the variants from the CDS regions
        self.vcf = self._get_variants_in_cds(vcf_tabix)
    
    def _get_variants_in_cds(self, vcf_tabix):
        """ gets the 1000 Genomes variants in the CDS regions of a gene
        """
        
        chrom = self.gene.get_chrom()
        
        vcf_records = []
        for start, end in self.gene.cds:
            if start != self.gene.cds_min:
                start -= 8
            if end != self.gene.cds_max:
                end += 8
            temp_records = vcf_tabix.fetch(str(chrom), start, end)
            
            # drop out any variants that already exist in vcf_records (mainly
            # affects CNVs and indels, which can be picked up multiple times if
            # they span more than one exon)
            duplicates_removed = []
            for var in temp_records:
                # check each var from the current exon against all the
                # previously added variants
                not_dup = True
                for previous in vcf_records:
                    if var.id == previous.id:
                        not_dup = False
                
                if not_dup:
                    duplicates_removed.append(var)
            
            vcf_records += duplicates_removed
        
        return vcf_records
    
    
