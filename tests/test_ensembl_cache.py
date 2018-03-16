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
    
    def test_retrieve_data(self):
        """ test that retrieve_data() works correctly
        """
        
        del self.cache.data
        # check that we raise an error if the data hasn't been set already
        with self.assertRaises(AttributeError):
            self.cache.retrieve_data()
        
        self.cache.data = "aaa"
        self.assertEqual(self.cache.retrieve_data(), "aaa")
    
    def test_get_key_from_url(self):
        """ test that get_key_from_url() works correctly
        """
        
        self.assertEqual(self.cache.get_key_from_url("http://beta.rest.ensembl.org/info/rest"), "info.rest")
        self.assertEqual(self.cache.get_key_from_url("http://beta.rest.ensembl.org/xrefs/symbol/homo_sapiens/ABO"), "xrefs.symbol.homo_sapiens.ABO")
        self.assertEqual(self.cache.get_key_from_url("http://beta.rest.ensembl.org/sequence/id/ENST00000378520?type=protein"), "sequence.id.ENST00000378520.protein")
        self.assertEqual(self.cache.get_key_from_url("http://beta.rest.ensembl.org/feature/id/ENSG00000175164?feature=transcript"), "feature.id.ENSG00000175164.transcript")
        self.assertEqual(self.cache.get_key_from_url("http://beta.rest.ensembl.org/sequence/id/ENST00000538324?type=genomic;expand_3prime=10;expand_5prime=10"), "sequence.id.ENST00000538324.genomic")
        self.assertEqual(self.cache.get_key_from_url("http://beta.rest.ensembl.org/sequence/id/ENST00000538324?type=cds"), "sequence.id.ENST00000538324.cds")
        self.assertEqual(self.cache.get_key_from_url("http://beta.rest.ensembl.org/feature/id/ENST00000538324?feature=exon"), "feature.id.ENST00000538324.exon")
        self.assertEqual(self.cache.get_key_from_url("http://beta.rest.ensembl.org/vep/human/id/rs3887873/consequences?"), "vep.human.id.rs3887873.consequences")
        self.assertEqual(self.cache.get_key_from_url("http://beta.rest.ensembl.org/vep/human/9:22125503-22125502:1/C/consequences?"), "vep.human.9_22125503-22125502_1.C.consequences")
    
    def test_check_if_data_in_cache(self):
        """ test that check_if_data_in_cache() works correctly
        """
        
        # set up the data to go in the database
        current_date = datetime.strftime(datetime.today(), "%Y-%m-%d")
        genome_build = "grch37"
        url = "http://beta.rest.ensembl.org/feature/id/temp1?feature=exon"
        key = self.cache.get_key_from_url(url)
        api_version = self.cache.api_version
        string = "temp_data"
        if IS_PYTHON3:
            string = string.encode("utf-8")
        temp_data = zlib.compress(string)
        
        if IS_PYTHON2:
            temp_data = buffer(temp_data)
        
        # check that the data not in the database returns False
        self.assertFalse(self.cache.check_if_data_in_cache(url))
        
        # insert the data in the database
        values = (key, genome_build, current_date, api_version, temp_data)
        self.cache.c.execute("INSERT into ensembl VALUES (?, ?, ?, ?, ?)", values)
        
        # check that the data in the database returns True
        self.assertTrue(self.cache.check_if_data_in_cache(url))
        
        # check that the cache's data object is set correctly if the row is in the database
        if IS_PYTHON2:
            self.assertEqual(self.cache.data, zlib.decompress(temp_data))
        elif IS_PYTHON3:
            self.assertEqual(self.cache.data, zlib.decompress(temp_data).decode("utf-8"))
        
        # check that obsolete data returns False
        old_date = datetime.strftime(datetime.today() - timedelta(days=181), "%Y-%m-%d")
        values = (key, genome_build, old_date, api_version, temp_data)
        self.cache.c.execute("REPLACE into ensembl VALUES (?, ?, ?, ?, ?)", values)
        self.assertFalse(self.cache.check_if_data_in_cache(url))
        
        # check that data from an old API version returns False
        old_api_version = "1.0.0"
        values = (key, genome_build, current_date, old_api_version, temp_data)
        self.cache.c.execute("REPLACE into ensembl VALUES (?, ?, ?, ?, ?)", values)
        self.assertFalse(self.cache.check_if_data_in_cache(url))
    
    def test_cache_url_data(self):
        """ test that cache_url_data works correctly
        """
        
        # set up the data to go in the database
        current_date = datetime.strftime(datetime.today(), "%Y-%m-%d")
        url = "http://beta.rest.ensembl.org/feature/id/temp2?feature=exon"
        key = self.cache.get_key_from_url(url)
        api_version = self.cache.api_version
        temp_data = "temp_data"
        
        # check that the data is not in before we insert it
        self.assertFalse(self.cache.check_if_data_in_cache(url))
        
        # insert the data, then check that it has gone in
        self.cache.cache_url_data(url, temp_data)
        self.assertTrue(self.cache.check_if_data_in_cache(url))
        
        # check that we can't cache the ensembl version string (so we have to
        # check this each time the script runs, rather than pulling from the
        # cache)
        url = "http://beta.rest.ensembl.org/info/rest"
        self.cache.cache_url_data(url, temp_data)
        self.assertFalse(self.cache.check_if_data_in_cache(url))
