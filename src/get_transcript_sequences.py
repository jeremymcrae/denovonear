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
        
        self.check_ensembl_api_version()
    
    def check_ensembl_api_version(self):
        """ check the ensembl api version matches a currently working version
        
        This function is included so when the api version changes, we notice the
        change, and we can manually check the responses for the new version.
        """
        
        headers = {"Content-Type": "application/json"}
        ext = "/info/rest"
        r = requests.get(self.server + ext, headers=headers)
        
        try:
            response = r.json()
            release = response["release"].split(".")
            
            major = release[0]
            minor = release[1]
            
            if major != "1" or minor != "6":
                raise ValueError("check ensembl api version")
            
            return
        except ValueError:
            self.rate_limit_ensembl_requests()
            self.check_ensembl_api_version()
    
    def ensembl_request(self, ext, sequence_id, headers):
        """ obtain sequence via the ensembl REST API
        """
        
        self.request_attempts += 1
        if self.request_attempts > 10:
            raise ValueError("too many attempts, figure out why its failing")
        
        self.rate_limit_ensembl_requests()
        try:
            r = requests.get(self.server + ext, headers=headers)
        except ConnectionError:
            self.ensembl_request(ext, sequence_id, headers)
        
        logging.warning("{0}\t{1}\t{2}\t{3}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), r.status_code, sequence_id, r.url))
        
        # we might end up passing too many requests per hour, just wait until
        # the period is finished before retrying
        if r.status_code == 429:
            time.sleep(int(r.headers["X-RateLimit-Reset"]))
            self.ensembl_request(ext, sequence_id, headers)
        # retry after 30 seconds if we get service unavailable error
        elif r.status_code == 503 or r.status_code == 504:
            time.sleep(30)
            self.ensembl_request(ext, sequence_id, headers)
        elif r.status_code != 200:
            raise ValueError("Invalid Ensembl response: " + str(r.status_code)\
                + " for " + str(sequence_id) + ". Submitted URL was: " + r.url)
        
        return r
    
    def get_genes_for_hgnc_id(self, hgnc_id):
        """ obatin the sensembl gene IDs that correspond to a HGNC symbol
        """
        
        headers = {"Content-Type": "application/json"}
        
        # http://beta.rest.ensembl.org/xrefs/symbol/homo_sapiens/KMT2A?content-type=application/json
        
        self.request_attempts = 0
        ext = "/xrefs/symbol/homo_sapiens/{0}".format(hgnc_id)
        r = self.ensembl_request(ext, hgnc_id, headers)
        
        genes = []
        for item in r.json():
            if item["type"] == "gene":
                genes.append(item["id"])
        
        return genes
        
    def get_transcript_ids_for_ensembl_gene_ids(self, gene_ids, hgnc_id):
    
        chroms = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", \
             "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", \
              "X", "Y"}
        
        headers = {"Content-Type": "application/json"}
        
        transcript_ids = []
        for gene_id in gene_ids:
            self.request_attempts = 0
            ext = "/feature/id/{0}?feature=transcript".format(gene_id)
            r = self.ensembl_request(ext, gene_id, headers)
            
            for item in r.json():
                # ignore transcripts not on the standard chromosomes 
                # (non-default chroms fail to map the known de novo variants 
                # to the gene location
                if item["Parent"] != gene_id or item["seq_region_name"] not in \
                        chroms or hgnc_id not in item["external_name"]:
                    continue
                transcript_ids.append(item["ID"])
        
        return transcript_ids
    
    def get_genomic_seq_for_transcript(self, transcript_id, expand):
        """ obtain the sequence for a transcript from ensembl
        """
        
        headers = {"Content-Type": "application/json"}
        
        self.request_attempts = 0
        ext = "/sequence/id/{0}?type=genomic;expand_3prime={1};expand_5prime={1}".format(transcript_id, expand)
        r = self.ensembl_request(ext, transcript_id, headers)
        
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
        r =  self.ensembl_request(ext, transcript_id, headers)
        
        return r.text
    
    def get_protein_seq_for_transcript(self, transcript_id):
        """ obtain the sequence for a transcript from ensembl
        """
        
        headers = {"Content-Type": "text/plain"}
        
        self.request_attempts = 0
        ext = "/sequence/id/" + transcript_id + "?type=protein"
        r =  self.ensembl_request(ext, transcript_id, headers)
        
        return r.text
    
    def get_chrom_for_transcript(self, transcript_id, hgnc_id):
        """ obtain the sequence for a transcript from ensembl
        """
        
        headers = {"Content-Type": "application/json"}
        
        self.request_attempts = 0
        ext = "/feature/id/" + transcript_id + "?feature=gene"
        r =  self.ensembl_request(ext, transcript_id, headers)
        
        for gene in r.json():
            if gene["external_name"] == hgnc_id:
                return gene["seq_region_name"]
        
        return None
    
    def get_exon_ranges_for_transcript(self, transcript_id):
        """ obtain the sequence for a transcript from ensembl
        """
        
        headers={"Content-Type": "application/json"}
        
        self.request_attempts = 0
        ext = "/feature/id/" + transcript_id + "?feature=exon"
        r = self.ensembl_request(ext, transcript_id, headers)
        
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
        r = self.ensembl_request(ext, transcript_id, headers)
        
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



    