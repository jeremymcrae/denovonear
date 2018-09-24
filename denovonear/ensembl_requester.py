""" Obtains genomic sequence, transcript IDs, gene, exon and CDS coordinates for
transcripts from Ensembl.
"""

from __future__ import print_function

import sys
import json
import time
import logging

IS_PYTHON3 = sys.version_info[0] == 3

if IS_PYTHON3:
    import urllib.request as request
    from urllib.error import HTTPError, URLError
else:
    import urllib2 as request
    from urllib2 import HTTPError, URLError

from denovonear.ensembl_cache import EnsemblCache

logging.basicConfig(filename='ensembl_requests.log', level=logging.WARNING)

class EnsemblRequest(object):
    """ Uses the Ensembl REST API to obtain gene information from Ensembl.
    
    Can find:
         - gene IDs for a HGNC symbol
         - transcript IDs for a gene ID
         - exon coordinates for an ensembl transcript ID
         - CDS coordinates for an ensembl transcript ID
         - transcript and genomic DNA sequences for an ensembl transcript ID
    """
    
    def __init__(self, cache_folder, genome_build):
        """ obtain the sequence for a transcript from ensembl
        
        Args:
            cache_folder: path to folder for caching data requested from Ensembl
            genome_build: string indicating the genome build ("grch37" or "grch38")
        """
        
        self.cache = EnsemblCache(cache_folder, genome_build)
        
        self.prior_time = time.time() - 1
        self.rate_limit = 0.067
        
        server_dict = {"grch37": "grch37.", "grch38": ""}
        
        self.server = "http://{}rest.ensembl.org".format(server_dict[genome_build])
        
        self.check_ensembl_api_version()
    
    def check_ensembl_api_version(self):
        """ check the ensembl api version matches a currently working version
        
        This function is included so when the api version changes, we notice the
        change, and we can manually check the responses for the new version.
        """
        
        self.attempt = 0
        headers = {"content-type": "application/json"}
        ext = "/info/rest"
        r = self.ensembl_request(ext, headers)
        
        response = json.loads(r)
        self.cache.set_ensembl_api_version(response["release"])
    
    def open_url(self, url, headers):
        """ open url with python libraries
        """
        
        data = self.cache.get_cached_data(url)
        if data is not None:
            return data, 200, headers
        
        self.rate_limit_ensembl_requests()
        req = request.Request(url, headers=headers)
        
        try:
            handler = request.urlopen(req)
        except HTTPError as error:
            # if we get a http error, we still process the status code, since a
            # later step deals with different status codes differently.
            handler = error
        except (URLError, ConnectionResetError, TimeoutError):
            # if we get a ConnectionResetError, assume something has gone wrong
            # with the server. Later code will wait before retrying.
            return '', 500, headers
        
        status_code = handler.getcode()
        response = handler.read()
        if IS_PYTHON3:
            response = response.decode("utf-8")
        
        # parse the headers into a key, value dictionary
        headers = dict(zip(map(str.lower, handler.headers.keys()), handler.headers.values()))
        
        now = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        logging.warning("{}\t{}\t{}".format(now, status_code, url))
        
        return response, status_code, headers
    
    def ensembl_request(self, ext, headers):
        """ obtain sequence via the ensembl REST API
        """
        
        self.attempt += 1
        if self.attempt > 5:
            raise ValueError("too many attempts, figure out why its failing")
        
        response, status, requested_headers = self.open_url(self.server + ext, headers=headers)
        
        # we might end up passing too many simultaneous requests, or too many
        # requests per hour, just wait until the period is finished before
        # retrying
        if status == 429:
            if "retry-after" in requested_headers:
                time.sleep(float(requested_headers["retry-after"]))
            elif "x-ratelimit-reset" in requested_headers:
                time.sleep(int(requested_headers["x-ratelimit-reset"]))
            
            return self.ensembl_request(ext, headers)
        # retry after 30 seconds if we get service unavailable error
        elif status in [500, 503, 504]:
            time.sleep(30)
            return self.ensembl_request(ext, headers)
        elif status != 200:
            raise ValueError("Invalid Ensembl response for {}\nheaders: {}\nresponse: {}".format(\
                    self.server + ext, requested_headers, response))
        
        # sometimes ensembl returns odd data. I don't know what it is, but the
        # json interpreter can't handle it. Rather than trying to catch it,
        # simply re-request the data
        if requested_headers["content-type"] == "application/json":
            try:
                json.loads(response)
            except ValueError:
                now = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
                logging.warning("{}\t{}\t{}\t{}\t{}".format(now,
                    status, self.server + ext,
                    "cannot obtain json output"))
                return self.ensembl_request(ext, requested_headers)
        
        self.cache.cache_url_data(self.server + ext, response)
        
        return response
    
    def get_genes_for_hgnc_id(self, hgnc_symbol):
        """ obtain the ensembl gene IDs that correspond to a HGNC symbol
        """
        
        headers = {"content-type": "application/json"}
        
        # http://grch37.rest.ensembl.org/xrefs/symbol/homo_sapiens/KMT2A?content-type=application/json
        
        self.attempt = 0
        ext = "/xrefs/symbol/homo_sapiens/{}".format(hgnc_symbol)
        r = self.ensembl_request(ext, headers)
        
        genes = []
        for item in json.loads(r):
            if item["type"] == "gene":
                genes.append(item["id"])
        
        return genes
    
    def get_previous_symbol(self, hgnc_symbol):
        """ sometimes we get HGNC symbols that do not match the ensembl rest version
        that we are currently using. We can look for earlier HGNC symbols for
        the gene using the service at rest.genenames.org
        
        Args:
            hgnc_symbol: HGNC symbol for the gene (eg "MLL2")
        
        Returns:
            list of deprecated gene symbols (eg ["KMT2A"])
        """
        
        ensembl_server = self.server
        gene_names_server = "http://rest.genenames.org"
        
        self.server = gene_names_server
        headers = {"accept": "application/json", "content-type": "application/json"}
        ext = "/fetch/symbol/{}".format(hgnc_symbol)
        try:
            r = self.ensembl_request(ext, headers)
        finally:
            self.server = ensembl_server
        
        gene_json = json.loads(r)
        
        prev_gene = []
        docs = gene_json["response"]["docs"]
        
        # strip out any gene entries that have been invalidated
        docs = [ x for x in docs if x["status"] != "Entry Withdrawn"]
        
        if len(docs) == 0:
            pass
        elif len(docs) > 1:
            raise ValueError("{0} has more than one alternate symbol, which I haven't accounted for.".format(hgnc_symbol))
        elif "prev_symbol" in docs[0]:
            prev_gene = docs[0]["prev_symbol"]
        
        return prev_gene
        
    def get_transcript_ids_for_ensembl_gene_ids(self, gene_ids, hgnc_symbols):
        """ fetch the ensembl transcript IDs for a given ensembl gene ID.
        
        Args:
            gene_ids: list of Ensembl gene IDs for the gene
            hgnc_symbols: list of possible HGNC symbols for gene
        """
        
        chroms = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", \
             "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", \
              "X", "Y"}
        
        headers = {"content-type": "application/json"}
        
        transcript_ids = []
        for gene_id in gene_ids:
            self.attempt = 0
            ext = "/overlap/id/{}?feature=transcript".format(gene_id)
            r = self.ensembl_request(ext, headers)
            
            for item in json.loads(r):
                # ignore non-coding transcripts
                if item["biotype"] not in ["protein_coding", "polymorphic_pseudogene"]:
                    continue
                
                # ignore transcripts not on the standard chromosomes
                # (non-default chroms fail to map the known de novo variants
                # to the gene location
                if item["Parent"] != gene_id or item["seq_region_name"] not in \
                        chroms or \
                        all([symbol not in item["external_name"] for symbol in hgnc_symbols]):
                    continue
                transcript_ids.append(item["id"])
        
        return transcript_ids
    
    def get_genomic_seq_for_transcript(self, transcript_id, expand):
        """ obtain the sequence for a transcript from ensembl
        """
        
        headers = {"content-type": "application/json"}
        
        self.attempt = 0
        ext = "/sequence/id/{0}?type=genomic;expand_3prime={1};expand_5prime={1}".format(transcript_id, expand)
        r = self.ensembl_request(ext, headers)
        
        gene = json.loads(r)
        
        seq = gene["seq"]
        seq_id = gene["id"]
        
        if seq_id != transcript_id:
            raise ValueError("ensembl gave the wrong transcript")
        
        desc = gene["desc"].split(":")
        chrom = desc[2]
        start = int(desc[3]) + expand
        end = int(desc[4]) - expand
        strand_temp = int(desc[5])
        
        strand = "+"
        if strand_temp == -1:
            strand = "-"
        
        return (chrom, start, end, strand, seq)
    
    def get_cds_seq_for_transcript(self, transcript_id):
        """ obtain the sequence for a transcript from ensembl
        """
        
        headers = {"content-type": "text/plain"}
        
        self.attempt = 0
        ext = "/sequence/id/{}?type=cds".format(transcript_id)
        
        return self.ensembl_request(ext,  headers)
    
    def get_protein_seq_for_transcript(self, transcript_id):
        """ obtain the sequence for a transcript from ensembl
        """
        
        headers = {"content-type": "text/plain"}
        
        self.attempt = 0
        ext = "/sequence/id/{}?type=protein".format(transcript_id)
        
        return self.ensembl_request(ext, headers)
    
    def get_genomic_seq_for_region(self, chrom, start_pos, end_pos):
        """ obtain the sequence for a genomic region
        """
        
        headers = {"content-type": "text/plain"}
        
        self.attempt = 0
        ext = "/sequence/region/human/{}:{}..{}:1".format(chrom, start_pos, end_pos)
        
        return self.ensembl_request(ext, headers)
    
    def get_chrom_for_transcript(self, transcript_id, hgnc_id):
        """ obtain the sequence for a transcript from ensembl
        """
        
        headers = {"content-type": "application/json"}
        
        self.attempt = 0
        ext = "/overlap/id/{}?feature=gene".format(transcript_id)
        r =  self.ensembl_request(ext, headers)
        
        for gene in json.loads(r):
            if gene["external_name"] == hgnc_id:
                return gene["seq_region_name"]
        
        return None
    
    def get_exon_ranges_for_transcript(self, transcript_id):
        """ obtain the sequence for a transcript from ensembl
        """
        
        headers = {"content-type": "application/json"}
        
        self.attempt = 0
        ext = "/overlap/id/{}?feature=exon".format(transcript_id)
        r = self.ensembl_request(ext, headers)
        
        exon_ranges = []
        for exon in json.loads(r):
            if exon["Parent"] != transcript_id:
                continue
            
            start = exon["start"]
            end = exon["end"]
            
            exon_ranges.append((start, end))
        
        return exon_ranges
    
    def get_cds_ranges_for_transcript(self, transcript_id):
        """ obtain the sequence for a transcript from ensembl
        """
        
        headers = {"content-type": "application/json"}
        
        self.attempt = 0
        ext = "/overlap/id/{}?feature=cds".format(transcript_id)
        r = self.ensembl_request(ext, headers)
        
        cds_ranges = []
        for cds_range in json.loads(r):
            if cds_range["Parent"] != transcript_id:
                continue
            
            start = cds_range["start"]
            end = cds_range["end"]
            
            cds_ranges.append((start, end))
        
        return cds_ranges
    
    def rate_limit_ensembl_requests(self):
        """ limit ensembl requests to one per 0.067 s
        """
        
        current_time = time.time()
        diff_time = current_time - self.prior_time
        
        # sleep until the current rate limit period is finished
        if diff_time < self.rate_limit:
            time.sleep(self.rate_limit - diff_time)
        
        # set the access time to now, for later checks
        self.prior_time = time.time()
