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
    
    def request_sequence(self, ext, sequence_id, headers):
        """ obtain sequence via the ensembl REST API
        """
        
        self.request_attempts += 1
        if self.request_attempts > 10:
            raise ValueError("too many attempts, figure out why its failing")
        
        self.rate_limit_ensembl_requests()
        r = requests.get(self.server + ext, headers=headers)
        
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
    
    def get_genomic_seq_for_transcript(self, transcript_id, expand):
        """ obtain the sequence for a transcript from ensembl
        """
        
        headers = {"Content-Type": "application/json"}
        
        self.request_attempts = 0
        ext = "/sequence/id/{0}?type=genomic;expand_3prime={1};expand_5prime={1}".format(transcript_id, expand)
        r = self.request_sequence(ext, transcript_id, headers)
        
        gene = r.json()
        
        seq = gene["seq"]
        seq_id = gene["id"]
        
        if seq_id != transcript_id:
            raise ValueError("ensembl gave the wrong transcript")
        
        desc = gene["desc"].split(":")
        chrom = desc[2]
        start = int(desc[3]) + expand
        end = int(desc[4]) - expand
        strand = desc[5]
        
        if int(strand) == -1:
            strand = "-"
        else:
            strand = "+"
        
        return (chrom, start, end, strand, seq)
    
    def get_cds_seq_for_transcript(self, transcript_id):
        """ obtain the sequence for a transcript from ensembl
        """
        
        headers = {"Content-Type": "text/plain"}
        
        self.request_attempts = 0
        ext = "/sequence/id/" + transcript_id + "?type=cds"
        r =  self.request_sequence(ext, transcript_id, headers)
        
        return r.text
    
    def get_protein_seq_for_transcript(self, transcript_id):
        """ obtain the sequence for a transcript from ensembl
        """
        
        headers = {"Content-Type": "text/plain"}
        
        self.request_attempts = 0
        ext = "/sequence/id/" + transcript_id + "?type=protein"
        r =  self.request_sequence(ext, transcript_id, headers)
        
        return r.text
    
    def get_exon_ranges_for_transcript(self, transcript_id):
        """ obtain the sequence for a transcript from ensembl
        """
        
        headers={"Content-Type": "application/json"}
        
        self.request_attempts = 0
        ext = "/feature/id/" + transcript_id + "?feature=exon"
        r = self.request_sequence(ext, transcript_id, headers)
        
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
        
        headers={"Content-Type": "application/json"}
        
        self.request_attempts = 0
        ext = "/feature/id/" + transcript_id + "?feature=cds"
        r = self.request_sequence(ext, transcript_id, headers)
        
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



    