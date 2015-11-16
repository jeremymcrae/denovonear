
# try interpro for domain information

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import json
import time
import logging
import platform

IS_PYTHON2 = sys.version_info[0] == 2
IS_PYTHON3 = sys.version_info[0] == 3

# load python version specific url opener library
if IS_PYTHON2:
    import urllib2 as request
    from urllib import urlencode
elif IS_PYTHON3:
    import urllib.request as request
    from urllib.parse import urlencode
else:
    raise ValueError("unknown python version")

from denovonear.gene_plot.interpro_parser import InterProParser

try:
    import xml.etree.cElementTree as et
except ImportError:
    import xml.etree.ElementTree as et

logging.basicConfig(filename='interpro_requests.log',level=logging.INFO)

class InterproRequest(InterProParser):
    """ Uses the InterPro REST API to obtain domain information from InterPro.
    """
    
    rate_limit = 0.067
    job_limit = 30
    title = "DomainPlotter"
    
    def __init__(self, email_address):
        """ obtain the sequence for a transcript from ensembl
        """
        
        self.prior_time = time.time() - 1
        self.submitted_jobs = 0
        self.results = []
        
        self.email = email_address
        self.server = "http://www.ebi.ac.uk/Tools/services/rest/iprscan5"
    
    def getUserAgent(self):
        """ User-agent for request (see RFC2616).
        """
        
        # Agent string for urllib2 library.
        urllib_agent = "Python-urllib {0}".format(request.__version__)
        clientVersion = '0'
        
        # Prepend client specific agent string.
        user_agent = "Interpro REST API ver{0}: Python {1} on {2} with {3}".format(
            clientVersion,
            platform.python_version(), platform.system(),
            urllib_agent
        )
        
        return user_agent
        
    def open_url(self, url, headers, post_data=False):
        """ open url with python libraries
        """
        
        req = request.Request(url, headers=headers)
        
        try:
            if post_data == False:
                handler = request.urlopen(req)
            else:
                data = urlencode(post_data)
                if IS_PYTHON3:
                    data = data.encode("utf-8")
                handler = request.urlopen(req, data)
        except Exception as e:
            handler = e
        
        status_code = handler.getcode()
        response = handler.read()
        if IS_PYTHON3:
            response = response.decode("utf-8")
        
        # parse the headers into a key, value dictionary
        headers = {}
        for key, value in zip(handler.headers.keys(), handler.headers.values()):
            headers[key.lower()] = value
        
        return response, status_code, headers
    
    def api_request(self, ext, headers, post_data=False):
        """ obtain sequence via the EBI InterPro REST API
        """
        
        self.request_attempts += 1
        if self.request_attempts > 5:
            raise ValueError("too many attempts, figure out why its failing")
        
        self.rate_limit_requests()
        response, status_code, requested_headers = self.open_url(self.server + ext, headers=headers, post_data=post_data)
        logging.info("{0}\t{1}\t{2}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), status_code, self.server + ext))
        
        if status_code != 200:
            raise ValueError("Invalid response: " + str(status_code)\
                + ".\nSubmitted URL was: " + self.server + ext + "\nheaders: " \
                + str(requested_headers) + "\nresponse: " + response)
        
        return response
    
    def rate_limit_requests(self):
        """ limit api requests to one per 0.335 s
        """
        
        sleeping = True
        while sleeping:
            if (time.time() - self.prior_time) > self.rate_limit:
                self.prior_time = time.time()
                sleeping = False
                
            time.sleep(0.01)
    
    def poll_until_job_finished(self, job_id):
        """ given a submitted job, periodically query until the job has finished
        """
        
        print("server is analysing: {0}".format(job_id))
        
        job_status = "PENDING"
        while job_status != "FINISHED":
            self.request_attempts = 0
            time.sleep(30)
            job_status = self.check_job_status(job_id)
            print(job_status)
            
            # quit out on various errors described at:
            # http://www.ebi.ac.uk/Tools/webservices/services/pfa/iprscan5_rest
            if job_status == "ERROR" or job_status == "FAILURE" or job_status == "NOT_FOUND":
                raise ValueError("submitted job has gone wrong: " + job_id)
            
        return self.get_result(job_id)
    
    def submit_interpro_run(self, protein_sequence):
        """ submit interpro run, to predict domains in a protein sequence
        """
        
        self.request_attempts = 0
        
        post_data = {}
        post_data["email"] = self.email
        post_data["title"] = self.title
        post_data["sequence"] = protein_sequence.rstrip() + "\n"
        
        header = {"User-Agent": self.getUserAgent()}
        
        job_id = self.api_request("/run/", header, post_data)
        self.job_id = job_id
        self.submitted_jobs += 1
        
        return self.poll_until_job_finished(job_id)
    
    def check_job_status(self, job_id):
        
        self.request_attempts = 0
        header = {"User-Agent": self.getUserAgent()}
        status = self.api_request("/status/{}".format(job_id), header)
        
        return status
    
    def get_result_types(self, job_id):
        
        self.request_attempts = 0
        header = {"User-Agent": self.getUserAgent()}
        result_types = self.api_request("/resulttypes/{}".format(job_id), header)
        
        root = et.fromstring(result_types)
        
        results = []
        for child in root:
            entry = {}
            entry["descriptor"] = child.find("description").text
            entry["identifier"] = child.find("identifier").text
            entry["file_suffix"] = child.find("fileSuffix").text
            entry["label"] = child.find("label").text
            entry["media_type"] = child.find("mediaType").text
            
            results.append(entry)
        
        return results
    
    def get_result(self, job_id):
        """ obtain the domain predictions for a InterPro REST API job ID
        """
        
        self.request_attempts = 0
        header = {"User-Agent": self.getUserAgent()}
        api_response = self.api_request("/result/{}/tsv".format(job_id), header)
        
        domains = self.parse_interpro_results(api_response)
        domains = self.standardise_names(domains)
        
        return domains
