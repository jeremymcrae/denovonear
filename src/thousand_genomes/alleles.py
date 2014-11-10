""" class to handle thousand genomes allele functions
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from src.thousand_genomes.sample_populations import sample_pops, populations

class Alleles(object):
    """ this class is intended to be a subclass of Extract1000Genomes, and
    contains functions to check allele frequencies in a 1000 Genomes VCF
    """
    
    _var_id = None
    geno_dict = {0: "hom_ref", 1: "het", 2: "hom_alt"}
    counts = None
    
    def _get_samples_with_allele(self, record, allele):
        """ get a list of sample IDs that have 
        
        Args:
            var: a VcfRecord for a variant
            allele: allele code that we wish to examine
        """
        
        match = set(["het", "hom_alt"])
        if allele == record.ref:
            match = set(["hom_ref", "het"])
        
        sample_n = 0
        sample_ids = set([])
        for data in record:
            split_data = data.split(":")
            
            if len(split_data) > 1 and split_data[0] != "GT":
                sample_id = record.samples[sample_n]
                population = sample_pops[sample_id]
                vcf_genotype = split_data[0]
                
                # females raise errors for variants on the Y-chromosome, just  
                # skip past those individuals, since they can't contribute to a 
                # genotype or population count
                try:
                    geno = self._convert_genotype(vcf_genotype)
                except ValueError:
                    continue
                
                sample_n += 1
                
                if geno in match:
                    sample_ids.add(sample_id)
        
        return (sample_ids, sample_n)
    
    def _get_major_and_minor_allele_by_freq(self, record):
        """ determines the major and minor allele from freq in 1000 genomes
        """
        
        self.counts = self._get_geno_counts(record)
        
        ref_count = 0
        alt_count = 0
        for population in self.counts:
            pop = self.counts[population]
            
            # count the reference and alternate alleles
            if "hom_ref" in pop:
                ref_count += 2 * pop["hom_ref"]
            
            if "het" in pop:
                ref_count += pop["het"]
                alt_count += pop["het"]
            
            if "hom_alt" in pop:
                alt_count += 2 * pop["hom_alt"]
            
        # get the alternate and reference allele frequencies
        try:
            ref_freq = ref_count/(ref_count + alt_count)
            alt_freq = alt_count/(ref_count + alt_count)
        except ZeroDivisionError:
            ref_freq = 0
            alt_freq = 0
        
        if ref_freq >= alt_freq:
            major_allele = record.ref
            minor_allele = record.alt[0]
        else:
            major_allele = record.alt[0]
            minor_allele = record.ref
        
        return major_allele, minor_allele
    
    def _get_max_maf(self, record):
        """ gets the highest 1000 genomes continental population minor allele 
        frequency for a single variant.
        
        Args:
            record: VcfRecord object
        
        Returns:
            the maximum minor allele frequency found amongst the 1000 genomes 
            continental populations.
        """
        
        self.counts = self._get_geno_counts(record)
        
        mafs = []
        for population in self.counts:
            pop = self.counts[population]
            ref_count = 0
            alt_count = 0
            
            # count the reference and alternate alleles
            if "hom_ref" in pop:
                ref_count += 2 * pop["hom_ref"]
            
            if "het" in pop:
                ref_count += pop["het"]
                alt_count += pop["het"]
            
            if "hom_alt" in pop:
                alt_count += 2 * pop["hom_alt"]
            
            maf = self._get_min_freq(ref_count, alt_count)
            mafs.append(maf)
        
        return max(mafs)
    
    def _get_max_minor_genotype_frequency(self, record):
        """ gets the highest 1000 genomes continental population minor  
        genotype frequency for a single variant.
        
        Args:
            record: VcfRecord object
        
        Returns:
            the maximum minor genotype frequency found amongst the 1000 genomes 
            continental populations.
        """
        
        self.counts = self._get_geno_counts(record)
        
        geno_freqs = []
        for population in self.counts:
            pop = self.counts[population]
            ref = 0
            alt = 0
            
            # count the reference and alternate alleles
            if "hom_ref" in pop:
                ref = pop["hom_ref"]
            
            if "hom_alt" in pop:
                alt = pop["hom_alt"]
            
            geno_freq = self._get_min_freq(ref, alt)
            geno_freqs.append(geno_freq)
        
        return max(geno_freqs)
    
    def _get_min_freq(self, ref_count, alt_count):
        """ get the minimum frequency from ref and alt counts
        
        This covers both allele and genotype frequencies, all we need is the
        count for the ref allele and the alt allele.
        
        Args:
            ref_count: int count for the reference
            alt_ccount: int count for the alternate
        
        Returns:
            the frequency (proportion) of the ref or alt (whichever is smaller)
        """
        
        # get the alternate and reference frequencies
        try:
            ref_freq = ref_count/(ref_count + alt_count)
        except ZeroDivisionError:
            ref_freq = 0
        
        alt_freq = 1 - ref_freq
        
        # the minor frequency will be the smaller of the alt and ref
        # frequencies (since the ref allele is merely the allele in the 
        # reference genome, not necessarily the most common allele in human 
        # populations)
        minor = min(ref_freq, alt_freq)
        
        return minor
    
    def _get_geno_counts(self, record):
        """ tally the different genotypes across the different populations
        
        Args:
            record: VcfRecord object
        
        Returns:
            dictionary of genotype counts by 1000 genomes continental population
        """
        
        # avoid checking the genotype counts more than once
        if record.id == self._var_id:
            return self.counts
        
        # avoid using the null "." ID, instead use the genome coordinate, which
        # is more likely to be unique
        self._var_id = record.id
        if record.id == ".":
            self._var_id = "{0}:{1}".format(record.contig, record.pos)
        
        counts = {}
        for population in populations:
            counts[population] = {"hom_ref": 0, "het": 0, "hom_alt": 0}
        
        pos = 0
        for data in record:
            split_data = data.split(":")
            
            if len(split_data) > 1 and split_data[0] != "GT":
                sample_id = record.samples[pos]
                population = sample_pops[sample_id]
                vcf_genotype = split_data[0]
                
                # females raise errors for variants on the Y-chromosome, just  
                # skip past those individuals, since they can't contribute to a 
                # genotype or population count
                try:
                    geno = self._convert_genotype(vcf_genotype)
                except ValueError:
                    continue
                
                counts[population][geno] += 1
                pos += 1
        
        return counts
        
    def _convert_genotype(self, vcf_genotype):
        """ converts a VcfRecord SNP genotype to genotype string
        
        Args:
            vcf_genotype: SNP genotype string e.g. "0|0", "1|0", or single 
                character genotype eg "0" or "1" for hemizygous individuals
        
        Returns:
            string value indicating the type of genotype ("hom_ref", "het", 
            "hom_alt)  
        """
        
        nonref_count = int(vcf_genotype[0]) + int(vcf_genotype[-1])
        
        return self.geno_dict[nonref_count]