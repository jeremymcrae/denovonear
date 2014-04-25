""" caches sequence information requested from ensembl, so we don't have to
repeatedly re-request information from the REST API if we have done so recently.
"""

import os
import json
from datetime import datetime

class EnsemblCache(object):
    """ Instead of repeatedly re-acquiring data from Ensembl each run, cache
    the requested data for faster retrieval
    """
    
    def __init__(self, cache_folder):
        """ initialise the class with the local cache folder
        
        Args:
            cache_folder: path to the cache 
        """
        
        self.cache_folder = cache_folder
    
    def set_ensembl_api_version(self, version):
        """ set the ensembl API version, so we can check for obsolete data
        
        Args:
            version: Ensembl API version string eg "1.6.0"
        """
        
        self.api_version = version
    
    def check_if_data_in_cache(self, url):
        """ checks if the data for a url has already been stored in the cache
        
        Makes sure that the data is not out of date, by only allowing data that 
        is less than six months old, and produced by the same Ensembl REST API 
        version.
        
        Args:
            url: URL for the Ensembl REST service
        
        Returns:
            True/False for having the data in the cache
        """
        
        path = self.parse_url(url)
        
        # if the data has been cached, check that it is not out of date, and 
        # the data was generated from the same Ensembl API version
        if os.path.exists(path):
            data = self.open_url_data(url)
            self.data = json.loads(data)
            api_version = self.data["ensembl_api"]
            
            cache_date = datetime.strptime(self.data["cache_date"], "%Y-%m-%d")
            diff = datetime.today() - cache_date
            
            if diff.days < 180 and self.api_version == api_version:
                return True
        
        return False
    
    def open_url_data(self, url):
        """ open the json data for a url, but don't do any processing
        
        Args:
            url: URL for the Ensembl REST service
        
        Returns:
            contents of a cached file for the URL.
        """
        
        path = self.parse_url(url)
        
        f = open(path, "r")
        data = f.read()
        f.close()
        
        return data
    
    def retrieve_data(self):
        """ retrieves data for a URL (checked earlier that this exists)
        
        The function is only used if the data exists in the cache. In order to 
        check if the data is in the cache, we loaded the file, and lightly 
        parsed the contents (to check if the data is obsolete). Rather than 
        reloading the file, just use the data loaded earlier.
        
        Returns:
            the data object for the URL (which will be parsed elsewhere)
        """
        
        return self.data["data"]
    
    def cache_url_data(self, url, data):
        """ cache the data retrieved from ensembl
        
        Args:
            url: URL for the Ensembl REST service
            data: response data from Ensembl
        """
        
        path = self.parse_url(url)
        
        # don't cache the ensembl version check
        if path.endswith("info/rest"):
            return
        
        # make sure a folder is availble to write a file into
        folder = os.path.dirname(path)
        self.make_folder(folder)
        
        current_date = datetime.strftime(datetime.today(), "%Y-%m-%d")
        
        # set up the data to write
        json_data = {}
        json_data["cache_date"] = current_date
        json_data["ensembl_api"] = self.api_version
        json_data["data"] = data
        
        # write the data to a file
        output = open(path, "w")
        json.dump(json_data, output)
        output.close()
    
    def make_folder(self, path):
        """ make a folder, starting from the first path that exists
        
        Args:
            path: folder path to save data into
        """
        
        folders = []
        
        # find the first folder in the path that exists
        while not os.path.exists(path):
            folders.insert(0, os.path.basename(path))
            path = os.path.dirname(path)
        
        # make each folder in turn, until we have completed the full path
        for folder in folders:
            path = os.path.join(path, folder)
            os.mkdir(path)
    
    def parse_url(self, url):
        """ parses the url into a list of folder locations
        
        Args:
            url: URL for the Ensembl REST service
        
        Returns:
            a parsed path to load from or save into
        """
        
        # Potential URLs types to parse into paths:
        # http://beta.rest.ensembl.org/info/rest
        # http://beta.rest.ensembl.org/xrefs/symbol/homo_sapiens/ABO
        # http://beta.rest.ensembl.org/sequence/id/ENST00000378520?type=protein
        # http://beta.rest.ensembl.org/feature/id/ENSG00000175164?feature=transcript
        # http://beta.rest.ensembl.org/sequence/id/ENST00000538324?type=genomic;expand_3prime=10;expand_5prime=10
        # http://beta.rest.ensembl.org/sequence/id/ENST00000538324?type=cds
        # http://beta.rest.ensembl.org/feature/id/ENST00000538324?feature=exon
        # http://beta.rest.ensembl.org/vep/human/id/rs3887873/consequences?
        # http://beta.rest.ensembl.org/vep/human/9:22125503-22125502:1/C/consequences?
        
        path = url.split("/")[3:]
        
        # fix the final bit of the url, which can have additional requirements,
        # none of which are necessary for uniquely defining the data
        path[-1] = path[-1].split(";")[0]
        
        # convert "feature=transcript" to "transcript" or "type=protein" to 
        # "protein"
        final = path[-1].split("?")
        if "=" in final[-1]:
            final[-1] = final[-1].split("=")[1]
        
        path = path[:-1] + final
        
        # replace characters not tolerated in paths
        for pos in range(len(path)):
            path[pos] = path[pos].replace(":", "_")
        
        path = os.path.join(self.cache_folder, *path)
        
        return path
            
        