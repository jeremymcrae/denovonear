"""
Copyright (c) 2015 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import unittest
import tempfile
import os
import sys
import shutil
import zlib
from datetime import datetime, timedelta
from threading import Thread, Event
import random
import hashlib
import time

from denovonear.ensembl_cache import EnsemblCache

IS_PYTHON2 = sys.version_info[0] == 2
IS_PYTHON3 = sys.version_info[0] == 3

class TestEnsemblCachePy(unittest.TestCase):
    """ unit test the EnsemblCache class
    """
    
    @classmethod
    def setUpClass(self):
        self.temp_dir = tempfile.mkdtemp()
        self.cache = EnsemblCache(self.temp_dir, "grch37")
        self.cache.set_ensembl_api_version("3.0.0")
    
    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.temp_dir)
    
    def test_get_key_from_url(self):
        """ test that get_key_from_url() works correctly
        """
        
        self.assertEqual(self.cache.get_key_from_url("http://rest.ensembl.org/info/rest"), "info.rest")
        self.assertEqual(self.cache.get_key_from_url("http://rest.ensembl.org/xrefs/symbol/homo_sapiens/ABO"), "xrefs.symbol.homo_sapiens.ABO")
        self.assertEqual(self.cache.get_key_from_url("http://rest.ensembl.org/sequence/id/ENST00000378520?type=protein"), "sequence.id.ENST00000378520.protein")
        self.assertEqual(self.cache.get_key_from_url("http://rest.ensembl.org/feature/id/ENSG00000175164?feature=transcript"), "feature.id.ENSG00000175164.transcript")
        self.assertEqual(self.cache.get_key_from_url("http://rest.ensembl.org/sequence/id/ENST00000538324?type=genomic;expand_3prime=10;expand_5prime=10"), "sequence.id.ENST00000538324.genomic")
        self.assertEqual(self.cache.get_key_from_url("http://rest.ensembl.org/sequence/id/ENST00000538324?type=cds"), "sequence.id.ENST00000538324.cds")
        self.assertEqual(self.cache.get_key_from_url("http://rest.ensembl.org/feature/id/ENST00000538324?feature=exon"), "feature.id.ENST00000538324.exon")
        self.assertEqual(self.cache.get_key_from_url("http://rest.ensembl.org/vep/human/id/rs3887873/consequences?"), "vep.human.id.rs3887873.consequences")
        self.assertEqual(self.cache.get_key_from_url("http://rest.ensembl.org/vep/human/9:22125503-22125502:1/C/consequences?"), "vep.human.9_22125503-22125502_1.C.consequences")
    
    def test_get_cached_data(self):
        """ test that get_cached_data() works correctly
        """
        
        # set up the data to go in the database
        url = "http://rest.ensembl.org/feature/id/temp1?feature=exon"
        string = "temp_data"
        
        # check that the data is not in the database to start
        self.assertIsNone(self.cache.get_cached_data(url))
        
        # insert the data in the database
        self.cache.cache_url_data(url, string)
        
        # check that some data is now in the database
        data = self.cache.get_cached_data(url)
        self.assertIsNotNone(data)
        
        # check that the data is correct if the row is in the database
        self.assertEqual(data, string)
    
    def test_get_cached_data_old_date(self):
        """ check that the cache ignores outdated data
        """
        url = "http://rest.ensembl.org/feature/id/temp1?feature=exon"
        string = "temp_data"
        today = datetime.today()
        long_ago = today - timedelta(days=181)
        
        # check that obsolete data returns False
        self.cache.today = long_ago
        self.cache.cache_url_data(url, string)
        self.assertIsNotNone(self.cache.get_cached_data(url))
        
        self.cache.today = today
        self.assertIsNone(self.cache.get_cached_data(url))
    
    def test_get_cached_data_old_api(self):
        """ check that data from an old API version returns False
        """
        url = "http://rest.ensembl.org/feature/id/temp1?feature=exon"
        string = "temp_data"
        old_version = "1.0.0"
        current_version = "2.0.0"
        
        self.cache.api_version = old_version
        self.cache.cache_url_data(url, string)
        
        # the data should be returns if it matches the current API version
        self.assertIsNotNone(self.cache.get_cached_data(url))
        
        # but no data if the API version for the data is outdated
        self.cache.api_version = current_version
        self.assertIsNone(self.cache.get_cached_data(url))
    
    def test_cache_url_data(self):
        """ test that cache_url_data works correctly
        """
        
        # set up the data to go in the database
        url = "http://rest.ensembl.org/feature/id/temp2?feature=exon"
        temp_data = "temp_data"
        
        # check that the data is not in before we insert it
        self.assertIsNone(self.cache.get_cached_data(url))
        
        # insert the data, then check that it has gone in
        self.cache.cache_url_data(url, temp_data)
        self.assertIsNotNone(self.cache.get_cached_data(url))
        
        # check that we can't cache the ensembl version string (so we check
        # this each time the script runs, rather than pulling from the cache)
        url = "http://rest.ensembl.org/info/rest"
        self.cache.cache_url_data(url, temp_data)
        self.assertIsNone(self.cache.get_cached_data(url))
    
    def test_cache_load(self):
        """ make sure the cache can handle a reasonable load
    
        This test uses multiple threads writing to the cache simultaneously to
        show the cache can handle the load. Failure is shown by an exception.
        """
    
        cache_dir = os.path.join(self.temp_dir, 'loading')
        os.mkdir(cache_dir)
        text = lambda l: '{:x}'.format(random.getrandbits(l * 4)).strip()
        url = lambda : 'example.com/base/sub/{}'.format(text(10))
        write = lambda cache: cache.cache_url_data(url(), text(100))
        
        class Runner(Thread):
            def __init__(self, counter=100):
                super(Runner, self).__init__()
                self.counter = counter
            def run(self):
                cache = EnsemblCache(cache_dir, 'grch37')
                cache.set_ensembl_api_version('6.0')
                while self.counter > 0:
                    write(cache)
                    self.counter -= 1
        
        try:
            threads = [ Runner() for x in range(50) ]
            [ x.start() for x in threads ]
            [ x.join() for x in threads ]
        except:
            self.fail("EnsemblCache failed under heavy load")
