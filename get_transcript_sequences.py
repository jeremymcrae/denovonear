""" Obtains genomic sequence for transcripts from Ensembl
"""

from __future__ import print_function
import time
import requests
import logging

logging.basicConfig(filename='ensembl_requests.log',level=logging.WARNING)

class GetTranscriptSequence(object):
    
    def __init__(self):
        """ obtain the sequence for a transcript from ensembl
        """
        
        self.prior_time = time.time()
        self.rate_limit = 0.335
        
        self.server = "http://beta.rest.ensembl.org"
        self.headers={"Content-Type": "text/plain"}
    
    def request_sequence(self, ext, sequence_id):
        """ obtain sequence via the ensembl REST API
        """
        
        self.request_attempts += 1
        if self.request_attempts > 10:
            raise ValueError("too many attempts, figure out why its failing")
        
        self.rate_limit_ensembl_requests()
        r = requests.get(self.server + ext, headers=self.headers)
        
        logging.warning("{0}\t{1}\t{2}\t{3}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), r.status_code, sequence_id, r.url))
        
        # we might end up passing too many requests per hour, just wait until
        # the period is finished before retrying
        if r.status_code == 429:
            time.sleep(int(r.headers["X-RateLimit-Reset"]))
            self.request_sequence(ext, sequence_id)
        # retry after 30 seconds if we get service unavailable error
        elif r.status_code == 503:
            time.sleep(30)
            self.request_sequence(ext, sequence_id)
        elif r.status_code != 200:
            raise ValueError("Invalid Ensembl response: " + str(r.status_code)\
                + " for " + sequence_id + ". Submitted URL was: " + r.url)
        
        return r
    
    def get_genomic_seq_for_transcript(self, transcript_id, expand_dist):
        """ obtain the sequence for a transcript from ensembl
        """
        
        self.request_attempts = 0
        ext = "/sequence/id/{0}?type=genomic;expand_3prime={1};expand_5prime={1}".format(transcript_id, expand_dist)
        # ext = "/sequence/id/" + transcript_id + "?type=genomic;expand_3prime=" + str(expand_dist) + ";expand_5prime=" + str(expand_dist)
        r = self.request_sequence(ext, transcript_id)
        
        return r.text
    
    def get_cds_seq_for_transcript(self, transcript_id):
        """ obtain the sequence for a transcript from ensembl
        """
        
        self.request_attempts = 0
        ext = "/sequence/id/" + transcript_id + "?type=cds"
        r =  self.request_sequence(ext, transcript_id)
        
        return r.text
    
    def get_protein_seq_for_transcript(self, transcript_id):
        """ obtain the sequence for a transcript from ensembl
        """
        
        self.request_attempts = 0
        ext = "/sequence/id/" + transcript_id + "?type=protein"
        r =  self.request_sequence(ext, transcript_id)
        
        return r.text
    
    def get_gene_start_and_end_for_transcript(self, transcript_id, gene_id):
        """ obtain the star and end points for a transcript from ensembl
        """
        
        self.headers={"Content-Type": "application/json"}
        
        self.request_attempts = 0
        ext = "/feature/id/" + transcript_id + "?feature=gene"
        r = self.request_sequence(ext, transcript_id)
        
        # an initial check, remove this once it works
        matching_genes = 0
        for gene in r.json():
            if gene["external_name"] == gene_id:
                matching_genes += 1
        
        if matching_genes > 1:
            raise ValueError("found too many gene matches")
        
        for gene in r.json():
            if gene["external_name"] == gene_id:
                start = gene["start"]
                end = gene["end"]
                strand = gene["strand"]
                chrom = gene["seq_region_name"]
                
                if strand == "-1":
                    strand = "-"
                else:
                    strand = "+"
                
                return (start, end, strand, chrom)
        
        raise ValueError("failed to find " + gene_id + " in ensembl request")
    
    def get_exon_ranges_for_transcript(self, transcript_id):
        """ obtain the sequence for a transcript from ensembl
        """
        
        self.headers={"Content-Type": "application/json"}
        
        self.request_attempts = 0
        ext = "/feature/id/" + transcript_id + "?feature=exon"
        r = self.request_sequence(ext, transcript_id)
        
        exon_ranges = []
        for exon in r.json():
            if exon["Parent"] != transcript_id:
                continue
            
            start = exon["start"]
            end = exon["end"]
            
            exon_ranges.append((start, end))
        
        return exon_ranges
    
    def get_cds_ranges_for_transcript(self, transcript_id):
        """ obtain the sequence for a transcript from ensembl
        """
        
        self.headers={"Content-Type": "application/json"}
        
        self.request_attempts = 0
        ext = "/feature/id/" + transcript_id + "?feature=cds"
        r = self.request_sequence(ext, transcript_id)
        
        cds_ranges = []
        for cds_range in r.json():
            if cds_range["Parent"] != transcript_id:
                continue
            
            start = cds_range["start"]
            end = cds_range["end"]
            
            cds_ranges.append((start, end))
        
        return cds_ranges
    
    def rate_limit_ensembl_requests(self):
        """ limit ensembl requests to one per 0.335 s
        """
        
        sleeping = True
        while sleeping:
            if (time.time() - self.prior_time) > self.rate_limit:
                self.prior_time = time.time()
                sleeping = False
                
            time.sleep(0.01)



    