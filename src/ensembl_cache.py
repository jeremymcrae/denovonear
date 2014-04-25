""" caches sequyence information requested from ensembl, so we don't have to
re-reqest information from the REST API if we have done so recently.
"""

from __future__ import print_function

import os
import json
from datetime import datetime

class EnsemblCache(object):
    """ Instead of repeatedly re-acquiring data from Ensembl each run, cache
    the requested data for faster retrieval
    """
    
    def __init__(self, cache_folder):
        """ initialise the class with the cache location
        """
        
        self.cache_folder = cache_folder
    
    def check_if_data_in_cache(self, url):
        """ checks if the data for a url has already been stored in the cache
        
        Makes sure that the data is not out of date, by only allowing data that 
        is less than six months old.
        """
        
        path = self.parse_url(url)
        
        # if the data has been cached, check that it is not out of date
        if os.path.exists(path):
            data = self.open_url_data(url)
            json_data = json.loads(data)
            
            cache_date = datetime.strptime(json_data["cache_date"], "%Y-%m-%d")
            diff = datetime.today() - cache_date
            
            if diff.days < 180:
                return True
        
        return False
    
    def open_url_data(self, url):
        """ open the json data for a url, but don't do any processing
        """
        
        path = self.parse_url(url)
        
        f = open(path, "r")
        data = f.readlines()
        f.close()
        
        return data
    
    def retrieve_url(self, url):
        """ retrieves data for a URL (checked earlier that this exists)
        """
        
        data = self.open_url_data(url)
        json_data = json.loads(data)
        
        return json_data["data"]
    
    def cache_url_data(self, url, data):
        """ save the data retrieved from ensembl
        """
        
        path = self.parse_url(url)
        current_date = datetime.strftime(datetime.today(), "%Y-%m-%d")
        
        json_data = {}
        json_data["cache_date"] = current_date
        json_data["data"] = data
        
        json.dump(json_data, path)
    
    def parse_url(self, url):
        """ parses the url into a list of folder locations
        """
        
        # http://beta.rest.ensembl.org/info/rest
        # http://beta.rest.ensembl.org/xrefs/symbol/homo_sapiens/ABO
        # http://beta.rest.ensembl.org/sequence/id/ENST00000378520?type=protein
        # http://beta.rest.ensembl.org/feature/id/ENSG00000175164?feature=transcript
        # http://beta.rest.ensembl.org/sequence/id/ENST00000538324?type=genomic;expand_3prime=10;expand_5prime=10
        # http://beta.rest.ensembl.org/sequence/id/ENST00000538324?type=cds
        # http://beta.rest.ensembl.org/feature/id/ENST00000538324?feature=exon
        # http://beta.rest.ensembl.org/vep/human/id/rs3887873/consequences?
        # http://beta.rest.ensembl.org/vep/human/9:22125503-22125502:1/C/consequences?
        
        path = url.lstrip("http://beta.rest.ensembl.org/")
        path = url.split("/")
        
        # fix the final bit of the url, which can have additional requirements,
        # none of which are necessary for uniquely defining the data
        path[-1] = path[-1].split(";")[0]
        
        # convert "feature=transcript" to "transcript" or "type=protein" to "protein"
        final = path[-1].split("?")
        if "=" in final[-1]:
            final[-1] = final[-1].split("=")[1]
        
        path = path[:-1] + final
        
        # replace characters not tolerated in paths
        for pos in range(len(path)):
            path[pos] = path[pos].replace(":", "_")
        
        path = os.path.join(self.cache_folder, path)
        
        return path
            
        