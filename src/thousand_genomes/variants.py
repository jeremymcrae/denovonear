""" script to obtain 1000 genomes variation data in CDS exons
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import pysam

from src.thousand_genomes.alleles import Alleles
from src.thousand_genomes.vcf import VCF
from src.thousand_genomes.functional_checks import FunctionalChecks
from src.site_specific_rates import SiteRates

class Extract1000Genomes(SiteRates, Alleles, VCF, FunctionalChecks):
    """ obtains data from 1000 Genomes VCF files
    """
    
    def __init__(self, vcf_folder="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/", frequency="allele"):
        """ intialise the class with the VCF location
        
        We can use VCFs on a ftp site, or a filesystem folder, since tabix copes
         with VCFs under either scenario.
        
        Args:
            vcf_folder: path or URL containing 1000 Genomes VCFs
            frequency: string, for whether to find minor allele, vs minor genotype frequencies
        """
            
        # set the url to the filesystem folder or ftp folder with the VCFs
        self.vcf_folder = vcf_folder
        
        if frequency == "allele":
            self._max_frequency = self._get_max_maf
        elif frequency == "genotype":
            self._max_frequency = self._get_max_minor_genotype_frequency
    
    def filter_variants(self, min_maf=0.01, max_maf=0.1, ignore_indels=True):
        """ filter for variants with suitable frequency, and functional properties
        """
        
        func = {"lof": [], "missense": [], "synonymous": []}
        for var in self.vcf:
            # ignore indels, CNVs and variants with > 1 alt allele
            if ignore_indels and (len(var.alt[0]) > 1 or len(var.ref) > 1):
                continue
            # ignore non-biallelic variants
            if len(var.alt[0]) > 1:
                continue
            
            # ignore extremely rare variants, or high-frequency variants.
            # TODO: also ignore singletons in 1000 genomes populations? i.e. 
            # TODO: a population might have one het individual, which takes the 
            # TODO: MAF between 0.01 and 0.1, which we would include, but
            if min_maf > self._max_frequency(var) or \
                    self._max_frequency(var) > max_maf:
                continue
            
            try:
                status = self.get_functional_status(var)
            except KeyError:
                continue
            
            (ref, alt) = self._get_major_and_minor_allele_by_freq(var)
            func[status].append((var, ref, alt))
        
        return func
    
    def find_samples_with_disrupted_alleles_in_gene(self):
        """ find the samples with disrupted allles in genes
        """
        
        samples = {"lof": set([]), "missense": set([])}
        sample_n = 0
        for var in self.vcf:
            # ignore non-biallelic variants
            if len(var.alt[0]) > 1:
                continue
            
            try:
                status = self.get_functional_status(var)
            except KeyError:
                continue
            
            if status == "synonymous":
                continue
            
            (ref, alt) = self._get_major_and_minor_allele_by_freq(var)
            sample_ids, temp_n = self.get_samples_with_allele(var, allele)
            samples[status] = samples[status] | sample_ids
            
            if temp_n > sample_n:
                sample_n = temp_n
        
        return samples, samples_n

