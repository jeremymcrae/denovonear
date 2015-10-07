""" caches sequence information requested from ensembl, so we don't have to
repeatedly re-request information from the REST API if we have done so recently.
"""

import os
import sqlite3
import sys
import time
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
        
        self.cache_folder = cache_folder
        self.cache_path = os.path.join(self.cache_folder, "ensembl_cache.db")
        
        self.genome_build = genome_build
        self.today = datetime.today()
        
        self.connect_to_database()
    
    def connect_to_database(self):
        """ connect to a sqlite database for the cached data (creates a new
        database if it doesn't already exist).
        """
        
        # create a folder for the cache database (this will raise an error
        # if we cannot create the folder)
        if not os.path.exists(self.cache_folder):
            os.mkdir(self.cache_folder)
        
        # generate a database with tables if it doesn't already exist
        if not os.path.exists(self.cache_path):
            conn = sqlite3.connect(self.cache_path)
            c = conn.cursor()
            try:
                c.execute(u"CREATE TABLE ensembl (key text PRIMARY KEY, genome_build text, cache_date text, api_version text, data blob)")
                conn.commit()
                c.close()
            except sqlite3.OperationalError:
                # occurs when multiple processes simultaneously try to create the
                # database
                
                # briefly sleep, so the process creating the database has time to
                # construct it
                c.close()
                time.sleep(5)
            
        self.conn = sqlite3.connect(self.cache_path)
        self.conn.row_factory = sqlite3.Row
        self.c = self.conn.cursor()
    
    def set_ensembl_api_version(self, version):
        """ set the ensembl API version, so we can check for obsolete data
        
        Args:
            version: Ensembl API version string eg "2.0.0"
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
        
        key = self.get_key_from_url(url)
        
        self.c.execute("SELECT * FROM ensembl WHERE key =? AND genome_build=?", (key, self.genome_build))
        row = self.c.fetchone()
        
        # if the data has been cached, check that it is not out of date, and
        # the data was generated from the same Ensembl API version
        if row is not None:
            api_version = row["api_version"]
            
            # by default extract the data, so we don't have to repeat the
            # sql query later.
            self.data = zlib.decompress(row["data"])
            if IS_PYTHON3:
                self.data = self.data.decode("utf-8")
            
            cache_date = datetime.strptime(row["cache_date"], "%Y-%m-%d")
            diff = self.today - cache_date
            
            if diff.days < 180 and self.api_version == api_version:
                return True
        
        return False
    
    def retrieve_data(self):
        """retrieves data for a URL (checked earlier that this exists)
        
        The function is only used if the data exists in the cache. In order to
        check if the data is in the cache, we selected database rows, and
        parsed the contents (to check if the data is obsolete). Rather than
        reselecting the data, just use the data loaded earlier.
        
        Returns:
            the data object for the URL (which will be parsed elsewhere)
        """
        
        return self.data
    
    def cache_url_data(self, url, data):
        """ cache the data retrieved from ensembl
        
        Args:
            url: URL for the Ensembl REST service
            data: response data from Ensembl
        """
        
        key = self.get_key_from_url(url)
        
        # don't cache the ensembl version check
        if key == "info.rest":
            return
        
        current_date = datetime.strftime(self.today, "%Y-%m-%d")
        
        # python3 zlib requires encoded strings
        if IS_PYTHON3:
            data = data.encode("utf-8")
        
        data = zlib.compress(data)
        
        # python2 sqlite3 can't write "8-bit bytestrings", but it can handle
        # buffer versions of the bytestrings
        if IS_PYTHON2:
            data = buffer(data)
        
        t = (key, self.genome_build, current_date, self.api_version, data)
        
        self.c.execute("INSERT OR REPLACE INTO ensembl (key, genome_build, cache_date, api_version, data) VALUES (?,?,?,?,?)", t)
        self.conn.commit()
    
    def get_key_from_url(self, url):
        """ parses the url into a list of folder locations
        
        Args:
            url: URL for the Ensembl REST service
        
        Returns:
            a parsed unique database key for the URLs data
        """
        
        key = url.split("/")[3:]
        
        # fix the final bit of the url, which can have additional requirements,
        # none of which are necessary for uniquely defining the data
        key[-1] = key[-1].split(";")[0]
        
        # convert "feature=transcript" to "transcript" or "type=protein" to
        # "protein"
        final = key[-1].split("?")
        if "=" in final[-1]:
            final[-1] = final[-1].split("=")[1]
        
        key = key[:-1] + final
        
        # replace characters not tolerated in keys
        for pos in range(len(key)):
            key[pos] = key[pos].replace(":", "_")
        
        # if the url ended with "?", then the final list element will be "",
        # which we should remove
        if key[-1] == "":
            key = key[:-1]
        
        key = ".".join(key)
        
        return key
