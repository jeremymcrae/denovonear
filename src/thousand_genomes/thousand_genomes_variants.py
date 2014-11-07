""" script to obtain 1000 genomes variation data in CDS exons
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os

import pysam

from src.thousand_genomes.thousand_genomes_sample_populations import sample_pops, populations
from src.site_specific_rates import SiteRates
from src.ensembl_requester import EnsemblRequest

class Extract1000Genomes(SiteRates):
    """ obtains data from 1000 Genomes VCF files
    """
    
    geno_dict = {0: "hom_ref", 1: "het", 2: "hom_alt"}
    
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
            self.__max_frequency = self.__get_max_maf
        elif frequency == "genotype":
            self.__max_frequency = self.__get_max_minor_genotype_frequency
    
    def set_gene(self, gene):
        """ set a gene for analysis - pull out CDS variants from the VCF
        
        Args:
            gene: Interval object, containing gene coordinates and sequence
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
        self.vcf = self.__get_variants_in_cds(vcf_tabix)
        
    def filter_variants(self, min_maf=0.01, max_maf=0.1, ignore_indels=True):
        """ filter for variants with suitable frequency, and functional properties
        """
        
        missense = []
        lof = []
        synonymous = []
        for var in self.vcf:
            # ignore indels, CNVs and variants with > 1 alt allele
            if ignore_indels and (len(var.alt[0]) > 1 or len(var.ref) > 1):
                continue
            # ignore non-biallelic variants
            if len(var.alt[0]) > 1:
                continue
            
            pos = var.pos
            if self.gene.strand == "-":
                pos -= 1
            
            self.counts = self.__get_geno_counts(var)
            # ignore extremely rare variants, or high-frequency variants.
            # TODO: also ignore singletons in 1000 genomes populations? i.e. 
            # TODO: a population might have one het individual, which takes the 
            # TODO: MAF between 0.01 and 0.1, which we would include, but
            if min_maf > self.__max_frequency(var) or \
                    self.__max_frequency(var) > max_maf:
                continue
            
            # get the distances to the closest exon boundaries
            exon_start, exon_end = self.gene.find_closest_exon(pos)
            self.boundary_dist = min(abs(exon_start - pos), abs(exon_end - pos))
            (ref, alt) = self.__get_major_and_minor_allele_by_freq(var)
            
            (lof_check, missense_check) = self.check_exon_spanning_deletion(pos, ref, exon_start, exon_end)
            if lof_check:
                lof.append((var, ref, alt))
                continue
            elif missense_check:
                missense.append((var, ref, alt))
                continue
            elif self.check_lof_indel(ref, alt):
                lof.append((var, ref, alt))
                continue
            elif self.check_missense_indel(ref, alt):
                missense.append((var, ref, alt))
                continue
            
            try:
                (lof_check, missense_check) = self.check_snv(pos, ref, alt)
            except KeyError:
                continue
            
            if lof_check:
                lof.append((var, ref, alt))
            elif missense_check:
                missense.append((var, ref, alt))
            else:
                synonymous.append((var, ref, alt))
        
        for (var, ref, alt) in lof:
            print(var.id)
        
        return (missense, lof, synonymous)
    
    def find_samples_with_lof_in_gene(self):
        
        lof_samples = set([])
        missense_samples = set([])
        sample_n = 0
        for var in self.vcf:
            # ignore non-biallelic variants
            if len(var.alt[0]) > 1:
                continue
            
            pos = var.pos
            if self.gene.strand == "-":
                pos -= 1
            
            # get the distances to the closest exon boundaries
            exon_start, exon_end = self.gene.find_closest_exon(pos)
            self.boundary_dist = min(abs(exon_start - pos), abs(exon_end - pos))
            (ref, alt) = self.__get_major_and_minor_allele_by_freq(var)
            
            (lof_check, missense_check) = self.check_exon_spanning_deletion(pos, ref, exon_start, exon_end)
            if lof_check:
                self.get_samples_with_allele(var, allele)
                continue
            elif missense_check:
                pass
                continue
            elif self.check_lof_indel(ref, alt):
                pass
                continue
            elif self.check_missense_indel(ref, alt):
                pass
                continue
            
            try:
                (lof_check, missense_check) = self.check_snv(pos, ref, alt)
            except KeyError:
                continue
            
            if lof_check:
                pass
            elif missense_check:
                pass
            else:
                pass
        
    
    def check_lof_indel(self, ref, alt):
        """ check if a variants is a frameshift indel
        """
        
        lof = False
        if len(alt) > 1 and (len(alt) - 1) % 3 != 0:
            lof = True
        if len(ref) > 1 and (len(ref) - 1) % 3 != 0:
            lof = True
        
        return lof
    
    def check_missense_indel(self, ref, alt):
        """ check if a variant is an in-frame indel
        """
        
        missense = False
        if len(alt) > 1 and (len(alt) - 1) % 3 == 0:
            missense = True
        if len(ref) > 1 and (len(ref) - 1) % 3 == 0:
            missense = True
        
        return missense
    
    def check_exon_spanning_deletion(self, pos, ref, exon_start, exon_end):
        """ check for deletions that span the exon boundaries
        
        I won't consider insertions, since insertions upstream of the exon
        probably wouldn't impact the transcript, whereas insertions that
        extend in distance past the exon end would simply extend the 
        transcript sequence, so it only matters if the insertion is in 
        frame or not. Obviously this will not always be true, eg
        insertions that introduce new intron/exon boundaries.
        """
         
        lof = False
        missense = False
        if pos < exon_start or pos + len(ref) > exon_end:
            ref_end = pos + len(ref)
            # check for deletions that surround an intron/exon boundary
            if pos < exon_start and ref_end > exon_start:
                lof = True
            elif pos < exon_end and ref_end > exon_end:
                lof = True
            # check for deletions that finish just short of an exon
            elif pos < exon_start and ref_end < exon_start:
                if abs(ref_end - exon_start) < 3:
                    lof = True
                elif abs(ref_end - exon_start) < 9:
                    missense = True
            # check for deletions that start just after an exon
            elif pos > exon_end:
                if self.boundary_dist < 3:
                    lof = False
                elif self.boundary_dist < 9:
                    missense = True
        
        return (lof, missense)
    
    def check_snv(self, pos, ref, alt):
        """ get the functional status of a SNV
        """
        
        # get the codon containing the variant, and the intra-codon position
        cds_pos = self.gene.convert_chr_pos_to_cds_positions(pos)
        codon_number = self.gene.get_codon_number_for_cds_position(cds_pos)
        codon_pos = self.gene.get_position_within_codon(cds_pos)
        codon = self.gene.get_codon_sequence(codon_number)
        
        # make sure the reference codon contains the reference allele (the
        # alleles might have been swapped since I select the more frequent 
        # allele as the reference, rather than directly from the reference
        # sequence)
        if codon[codon_pos] == alt:
            codon[codon_pos] = ref
        
        # figure out if the variant impacts the amino acid.
        # Occasionally the variant occurs in an incomplete transcript (eg 
        # rs28394186 in ENST00000336769). Those may give odd length codons
        # which we shall ignore (since they are so rare).
        initial_aa = self.gene.translate_codon(codon)
        mutated_aa = self.get_mutated_aa(alt, codon, codon_pos)
        
        if self.gene.strand == "-":
            pos += 1
        
        lof = False
        missense = False
        if self.missense_and_splice_region_check(initial_aa, mutated_aa, pos):
            missense = True
        elif self.loss_of_function_check(initial_aa, mutated_aa, pos):
            lof = True
        
        return (lof, missense)
    
    def get_samples_with_allele(self, record, allele):
        """ get a list of sample IDs that have 
        
        Args:
            var: a VcfRecord for a variant
            genotype: allele code that we wish to examine
        """
        
        match = set([1, 2])
        if allele == record.ref:
            match = set([0, 1])
        
        sample_n = 0
        sample_ids = []
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
                    geno = self.__convert_genotype(vcf_genotype)
                except ValueError:
                    continue
                
                sample_n += 1
                
                if geno in match:
                    sample_ids.append(sample_id)
        
        return (sample_ids, sample_n)
    
    def __get_variants_in_cds(self, vcf_tabix):
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
    
    def __get_major_and_minor_allele_by_freq(self, record):
        """ determines the major and minor allele from freq in 1000 genomes
        """
        
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
    
    def __get_max_maf(self, record):
        """ gets the highest 1000 genomes continental population minor allele 
        frequency for a single variant.
        
        Args:
            record: VcfRecord object
        
        Returns:
            the maximum minor allele frequency found amongst the 1000 genomes 
            continental populations.
        """
        
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
            
            maf = self.get_min_freq(ref_count, alt_count)
            mafs.append(maf)
        
        return max(mafs)
    
    def __get_max_minor_genotype_frequency(self, record):
        """ gets the highest 1000 genomes continental population minor  
        genotype frequency for a single variant.
        
        Args:
            record: VcfRecord object
        
        Returns:
            the maximum minor genotype frequency found amongst the 1000 genomes 
            continental populations.
        """
        
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
            
            geno_freq = self.get_min_freq(ref, alt)
            geno_freqs.append(geno_freq)
        
        return max(geno_freqs)
    
    def get_min_freq(self, ref_count, alt_count):
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
    
    def __get_geno_counts(self, record):
        """ tally the different genotypes across the different populations
        
        Args:
            record: VcfRecord object
        
        Returns:
            a dictionary of the different genotype counts by 1000 genomes 
            continental population
        """
        
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
                    geno = self.__convert_genotype(vcf_genotype)
                except ValueError:
                    continue
                
                counts[population][geno] += 1
                pos += 1
        
        return counts
        
    def __convert_genotype(self, vcf_genotype):
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

