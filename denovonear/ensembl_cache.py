""" caches sequence information requested from ensembl, so we don't have to
repeatedly re-request information from the REST API if we have done so recently.
"""

import os
import sqlite3
import sys
import time
import random
import zlib
from pathlib import Path
import json
from datetime import datetime
from urllib.parse import urlparse

class EnsemblCache(object):
    """ Instead of repeatedly re-acquiring data from Ensembl each run, cache
    the requested data for faster retrieval
    """
    today = datetime.today()
    def __init__(self, cache_folder):
        """ initialise the class with the local cache folder
        
        Args:
            cache_folder: path to the cache
        """
        cache_folder = Path(cache_folder)
        if not cache_folder.exists():
            cache_folder.mkdir()
        
        # generate a database with tables if it doesn't already exist
        path = cache_folder / "ensembl_cache.db"
        if not path.exists():
            try:
                with sqlite3.connect(str(path)) as conn:
                    with conn as cursor:
                        cursor.execute("CREATE TABLE ensembl " \
                            "(key text PRIMARY KEY, genome_build text, " \
                            "cache_date text, api_version text, data blob)")
            except sqlite3.OperationalError:
                time.sleep(random.uniform(1, 5))
        
        self.conn = sqlite3.connect(str(path))
        self.conn.row_factory = sqlite3.Row
        self.cursor = self.conn.cursor()
    
    def get_cached_data(self, url):
        """ get cached data for a url if stored in the cache and not outdated
        
        Args:
            url: URL for the Ensembl REST service
        
        Returns:
            data if data in cache, else None
        """
        key, build = self.get_key_from_url(url)
        self.cursor.execute("SELECT * FROM ensembl WHERE key=? AND genome_build=?",
            (key, build))
        row = self.cursor.fetchone()
        
        # if the data has been cached, check that it is not out of date
        if row is not None:
            data = zlib.decompress(row["data"])
            diff = self.today - datetime.strptime(row["cache_date"], "%Y-%m-%d")
            if diff.days < 180:
                return data
        
        return None
    
    def cache_url_data(self, url, data, attempt=0):
        """ cache the data retrieved from ensembl
        
        Args:
            url: URL for the Ensembl REST service
            data: response data from Ensembl, in bytes form
        """
        if attempt > 5:
            raise ValueError('too many attempts at writing to the cache')
        
        key, build = self.get_key_from_url(url)
        current_date = datetime.strftime(self.today, "%Y-%m-%d")
        
        compressed = zlib.compress(data)
        t = (key, build, current_date, '9', compressed)
        
        cmd = "INSERT OR REPLACE INTO ensembl " \
            "(key, genome_build, cache_date, api_version, data) VALUES (?,?,?,?,?)"
        try:
            with self.conn as conn:
                conn.execute(cmd, t)
        except sqlite3.OperationalError:
            # if we hit a sqlite locking error, wait a random time so conflicting
            # instances are less likely to reconflict, then retry
            time.sleep(random.uniform(1, 10))
            self.cache_url_data(url, data, attempt + 1)
    
    def get_key_from_url(self, url):
        """ parses the url into a list of folder locations
        
        We take a URL like:
            http://rest.ensembl.org/sequence/id/ENST00000538324?type=genomic;expand_3prime=10
        and turn it into 'sequence.id.ENST00000538324.genomic'
        
        Args:
            url: URL for the Ensembl REST service
        
        Returns:
            a parsed unique database key for the URLs data
        """
        url = urlparse(url)
        path = url.path.strip('/').replace('/', '.')
        
        # find the build from the url, but convert the default server to a build
        build = url.netloc.split('.')[0]
        build = 'grch38' if build == 'rest' else build
        
        # convert "LONG_ID?feature=transcript" to ['LONG_ID', "transcript"] etc
        suffix = url.query.split(';')[0]
        if "=" in suffix:
            _, suffix = suffix.split("=")
        
        key = path + '.' + suffix if suffix != '' else path
        # replace characters not tolerated in keys
        return key.replace(':', '_'), build
