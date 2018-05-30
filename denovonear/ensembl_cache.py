""" caches sequence information requested from ensembl, so we don't have to
repeatedly re-request information from the REST API if we have done so recently.
"""

import os
import sqlite3
import sys
import time
import random
import zlib
from datetime import datetime

IS_PYTHON2 = sys.version_info[0] == 2
IS_PYTHON3 = sys.version_info[0] == 3

class EnsemblCache(object):
    """ Instead of repeatedly re-acquiring data from Ensembl each run, cache
    the requested data for faster retrieval
    """
    
    def __init__(self, cache_folder, genome_build):
        """ initialise the class with the local cache folder
        
        Args:
            cache_folder: path to the cache
        """
        self.api_version = ('1')
        self.genome_build = genome_build
        self.today = datetime.today()
        
        if not os.path.exists(cache_folder):
            os.mkdir(cache_folder)
        
        # generate a database with tables if it doesn't already exist
        path = os.path.join(cache_folder, "ensembl_cache.db")
        if not os.path.exists(path):
            try:
                with sqlite3.connect(path) as conn:
                    with conn as cursor:
                        cursor.execute("CREATE TABLE ensembl " \
                            "(key text PRIMARY KEY, genome_build text, " \
                            "cache_date text, api_version text, data blob)")
            except sqlite3.OperationalError:
                time.sleep(random.uniform(1, 5))
        
        self.conn = sqlite3.connect(path)
        self.conn.row_factory = sqlite3.Row
    
    def set_ensembl_api_version(self, version):
        """ set the ensembl API version, so we can check for obsolete data
        
        Args:
            version: Ensembl API version string eg "2.0.0"
        """
        
        self.api_version = version
    
    def get_cached_data(self, url):
        """ get cached data for a url if stored in the cache and not outdated
        
        Args:
            url: URL for the Ensembl REST service
        
        Returns:
            data if data in cache, else None
        """
        
        key = self.get_key_from_url(url)
        
        with self.conn as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT * FROM ensembl WHERE key=? AND genome_build=?",
                (key, self.genome_build))
        row = cursor.fetchone()
        
        # if the data has been cached, check that it is not out of date, and
        # the data was generated from the same Ensembl API version
        if row is not None:
            api_version = row["api_version"]
            data = zlib.decompress(row["data"])
            if IS_PYTHON3:
                data = data.decode("utf-8")
            
            date = datetime.strptime(row["cache_date"], "%Y-%m-%d")
            diff = self.today - date
            
            if diff.days < 180 and self.api_version == api_version:
                return data
        
        return None
    
    def cache_url_data(self, url, data, attempt=0):
        """ cache the data retrieved from ensembl
        
        Args:
            url: URL for the Ensembl REST service
            data: response data from Ensembl
        """
        
        if attempt > 5:
            raise ValueError('too many attempts at writing to the cache')
        
        key = self.get_key_from_url(url)
        
        # don't cache the ensembl version check
        if key == "info.rest":
            return
        
        current_date = datetime.strftime(self.today, "%Y-%m-%d")
        
        # python3 zlib requires encoded strings
        if IS_PYTHON3:
            data = data.encode("utf-8")
        
        compressed = zlib.compress(data)
        
        # python2 sqlite3 can't write "8-bit bytestrings", but it can handle
        # buffer versions of the bytestrings
        if IS_PYTHON2:
            compressed = buffer(compressed)
        
        t = (key, self.genome_build, current_date, self.api_version, compressed)
        
        cmd = "INSERT OR REPLACE INTO ensembl " \
            "(key, genome_build, cache_date, api_version, data) VALUES (?,?,?,?,?)"
        try:
            with self.conn as cursor:
                cursor.execute(cmd, t)
        except sqlite3.OperationalError:
            # if we hit a sqlite locking error, wait a random time so conflicting
            # instances are less likely to reconflict, then retry
            time.sleep(random.uniform(1, 10))
            self.cache_url_data(url, data.decode('utf-8'), attempt + 1)
    
    def get_key_from_url(self, url):
        """ parses the url into a list of folder locations
        
        We take a URL like:
            http://rest.ensembl.org/sequence/id/ENST00000538324?type=genomic;expand_3prime=10;expand_5prime=10
        and turn it into 'sequence.id.ENST00000538324.genomic'
        
        Args:
            url: URL for the Ensembl REST service
        
        Returns:
            a parsed unique database key for the URLs data
        """
        
        key = url.split("/")[3:]
        
        # fix the final bit of the url, none of which uniquely define the data
        suffix = key.pop()
        suffix = suffix.split(";")[0]
        
        # convert "LONG_ID?feature=transcript" to ['LONG_ID', "transcript"] etc
        id = suffix.split("?", 1)
        suffix = id.pop()
        if "=" in suffix:
            _, suffix = suffix.split("=")
        
        key += id + [suffix]
        
        # replace characters not tolerated in keys and remove blank entries
        key = ( x.replace(':', '_') for x in key )
        key = ( x for x in key if x != '' )
        
        return ".".join(key)
