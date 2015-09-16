""" class to extract variant consequences from the Ensembl REST API
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from denovonear.ensembl_requester import EnsemblRequest

class EnsemblWithVariants(EnsemblRequest):
    """ class to extract variant consequences from the Ensembl REST API
    """
    
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
