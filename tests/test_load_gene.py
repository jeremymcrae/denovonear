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

import asyncio
import unittest
import tempfile
from pathlib import Path

from denovonear.load_gene import cds_length, construct_gene_object, \
    get_de_novos_in_transcript, get_transcript_ids, load_gene, \
    count_de_novos_per_transcript, minimise_transcripts
from denovonear.transcript import Transcript
from denovonear.rate_limiter import RateLimiter

async def call(func, *args, **kwargs):
    ''' call ensembl rest API function
    '''
    async with RateLimiter(15) as ensembl:
        return await func(ensembl, *args, **kwargs)

def _run(func, *args, **kwargs):
    return asyncio.get_event_loop().run_until_complete(call(func, *args, **kwargs))

class TestLoadGenePy(unittest.TestCase):
    """ unit test functions to load genes
    """
    
    def set_transcript(self):
        """ construct a transcript for a known gene
        """
        
        exon_ranges=[(120933859, 120934019), (120934219, 120934356),
            (120935876, 120936296)]
        cds_ranges=[(120934225, 120934356), (120935876, 120936013)]
        
        expected = Transcript("ENST00000242577", '12', 120933859, 120936296, "+")
        expected.set_exons(exon_ranges, cds_ranges)
        expected.set_cds(cds_ranges)
        
        cds = "ATGTGCGACCGAAAGGCCGTGATCAAAAATGCGGACATGTCGGAAGAGATGCAACAGGACTC" \
            "GGTGGAGTGCGCTACTCAGGCGCTGGAGAAATACAACATAGAGAAGGACATTGCGGCTCATATC" \
            "AAGAAGGAATTTGACAAGAAGTACAATCCCACCTGGCATTGCATCGTGGGGAGGAACTTCGGTA" \
            "GTTATGTGACACATGAAACCAAACACTTCATCTACTTCTACCTGGGCCAAGTGGCCATTCTTCT" \
            "GTTCAAATCTGGTTAA"
        
        genomic = "GGGCGGGGCCGGCGGGAGCAGGGCGGGGCCTGAGCACTAGGCGGCGGCGGCTGGCGTGGG" \
            "GCTGCTTAGATGCGCCACGGTTTCGGTAGCGACGGTATCTCTAGCCGGGCCTGAGCTGTGCTAGCA" \
            "CCTCCCCCAGGAGACCGTTGCAGTCGGCCAGCCCCCTTCTCCACGGTGAGAAACTCGGGGGGCCAG" \
            "GGGGTGTCCTCGCTGCCTTATTTCGCCCCACTCCGGACTTAGCCCTCCGCGTAGCCCGCGCTTCCT" \
            "GAGAAGTGGGGTGGGGGGCGTCGTCCCGTGGTGGCGCCGGCCGGGGTGGGGGCAGTTAGTGCCTGG" \
            "GGGGCGCGGCCCAACTCAACCCCTTACCCCAGGCCTTGCCCACTAGGTAACCATGTGCGACCGAAA" \
            "GGCCGTGATCAAAAATGCGGACATGTCGGAAGAGATGCAACAGGACTCGGTGGAGTGCGCTACTCA" \
            "GGCGCTGGAGAAATACAACATAGAGAAGGACATTGCGGCTCATATCAAGAAGGTGAGGATGGGCGC" \
            "GGGGGCCGATACGCAGCCGGGAGCAGGGGGTTCCTTCCCCCCGATCCTGCTTTCCTAAGGGCGCCT" \
            "GACAGGTCCCGGGAATACTGCTGGCGGCTTGGGGCGTAGAAGCTTCCAGAAAGGACGCAGATGCAT" \
            "TTTGCGCTCCTGTGGAGAAGACCAGACCCCCGGCGTCCGAAGTTTTTTTTTTTTTTTTTTTAATTA" \
            "CCCAGCTCCGCGGGGGGAAAGCGCCACCTAGCAACGGTATCTAAGATCAGGGAGCAGCGGTTCCCC" \
            "CTTCTGTGTGGTTCCTGCGCCGAGGATCCATCTGGGTGTTCCGGAGGGGGGAGCTGCGTGGGTGTT" \
            "TCCAGCCGGGCCGGGAGGAGATCTTGCCAGCCTTCCAGTGGGGAGTTGAGGGAAGGTGGTGGGTGG" \
            "TGGCGGGGCTGGGGGCTGGGGTAGGGGCTTGGTAAATGGCAGTCTAGAAAGCCGGCAGGACTGCCA" \
            "ACTTCTCGAGCAGTGTTTGCTGGAAGGGAAGAAAGCTGGCAGCCTAAGCCGTGGGAGGGTTCCAGT" \
            "CGAGAATGGGAAGATGAAAGACTTCAGATGGAACAGAAATAAATGCCTTTTTTGACAAACGCAGCA" \
            "GTGCGTGCCTCTAGCTTGCAAGAGCGTTACTCCCCTTCATAGCTTTAAAAGGTTTTCGCACTGCGT" \
            "GCAGTTAGAGTAGCTAAATCTTGTGTGACGCTCCACAAACACTTGTAAGAATTTTGCAGAGAAAGA" \
            "TAACCGTTGCCACCCAATGCCCCCCACAGGCATTCTACTCCCCAGTACCTCTTAGGGTGGGAGAAA" \
            "TGGTGAAGAGTTGTTCCTACAACTTGCTAACCTAGTGGACAGGGTAGTAGATTAGCATCATCCGGA" \
            "TAGATGTGAAGAGGACGGCTGTTTGGATAATAATTAAGGATAAAATTTGGCCAGTTGACAGATTCT" \
            "GTTTCCAGCAGTTTTTACAGCAACAGTGGAGTGCTTCAGTATTGTGTTCCTGTAAATTTAATTTTG" \
            "ATCCGCAATCATTTGGTATACAATGCTGTTTGAAGTTTTGTCCTATTGGAAAAGTCTTGTGTTGCA" \
            "GGGGTGCAGTTAAGATCTTTGTGATGAGGAATGGGATGGGCTAATTTTTTGCCGTTTTCTTGGAAT" \
            "TGGGGGCATGGCAAATACAGTAGGGTAGTTTAGTTCTCTACACAGAACATGATAAACTACACCTGT" \
            "TGATGTCACCGTCTGTCAATGAATATTATAGAAGGTATGAAGGTGTAATTACCATAATAACAAAAC" \
            "ACCCTGTCTTTAGGGCTGACCTTTCGTCCTTTGACCTCCTCAGCCTCCATTCCCATCTTCGCTCAG" \
            "ACTGCAAGTATGTTTGTATTAATGTACTATGTAGGCGGCTTGGAGCTGGGGAACATTCTTTCATTC" \
            "TAAGAATTTGCAGATGCTGACGTTCCTCCTTTCTGCCCCTACAGGCTCTGGCTTATCCAAGAGGCA" \
            "AACACTGACCTCTGGTAATTAAAATCCTAGTTCTTTTCTTTTGTCTTTTCCAGGAATTTGACAAGA" \
            "AGTACAATCCCACCTGGCATTGCATCGTGGGGAGGAACTTCGGTAGTTATGTGACACATGAAACCA" \
            "AACACTTCATCTACTTCTACCTGGGCCAAGTGGCCATTCTTCTGTTCAAATCTGGTTAAAAGCATG" \
            "GACTGTGCCACACACCCAGTGATCCATCCAAAAACAAGGACTGCAGCCTAAATTCCAAATACCAGA" \
            "GACTGAAATTTTCAGCCTTGCTAAGGGAACATCTCGATGTTTGAACCTTTGTTGTGTTTTGTACAG" \
            "GGCATTCTCTGTACTAGTTTGTCGTGGTTATAAAACAATTAGCAGAATAGCCTACATTTGTATTTA" \
            "TTTTCTATTCCATACTTCTGCCCACGTTGTTTTCTCTCAAAATCCATTCCTTTAAAAAATAAATCT" \
            "GATGCAGATGTGTATGTGTGTG"
        
        expected.add_cds_sequence(cds)
        expected.add_genomic_sequence(genomic, offset=10)
        
        return expected
    
    def _construct_transcript(self, tx_id):
        ''' helper function to prepare a Transcript given a transcript ID
        '''
        return _run(construct_gene_object, tx_id)
    
    def test_cds_length(self):
        ''' check we get the correct CDS length for a transcript
        '''
        tx = self.set_transcript()
        self.assertEqual(cds_length(tx), 270)
        
        tx = self._construct_transcript('ENST00000242577')
        self.assertEqual(cds_length(tx), 270)
        
        tx = self._construct_transcript('ENST00000392509')
        self.assertEqual(cds_length(tx), 270)
        
        tx = self._construct_transcript('ENST00000549649')
        self.assertEqual(cds_length(tx), 126)
    
    def test_construct_gene_object(self):
        """
        """
        
        transcript_id = "ENST00000242577"
        transcript = _run(construct_gene_object, transcript_id)
        
        expected = self.set_transcript()
        
        self.assertEqual(transcript, expected)
        self.assertEqual(transcript.get_genomic_sequence(), expected.get_genomic_sequence())
        self.assertEqual(transcript.get_cds_sequence(), expected.get_cds_sequence())
    
    def test_get_de_novos_in_transcript(self):
        """ test that we can identify de novos within the CDS of a transcript
        """
        
        exon_ranges = [(10, 20), (30, 40), (90, 100)]
        cds_ranges = [(30, 40), (90, 95)]
        
        # define a simple transcript
        tx = Transcript("test1", '1', 10, 100, "+")
        tx.set_exons(exon_ranges, cds_ranges)
        tx.set_cds(cds_ranges)
        
        # check that only the site in the CDS is returned
        sites = [15, 35, 100]
        self.assertEqual(get_de_novos_in_transcript(tx, sites), [35])
        
        # check that we can return multiple sites in the CDS
        sites = [15, 35, 90]
        self.assertEqual(get_de_novos_in_transcript(tx, sites), [35, 90])
        
        # check if we pass in an empty list, we get one back
        self.assertEqual(get_de_novos_in_transcript(tx, []), [])
    
    def test_get_transcript_ids(self):
        """ check that we get the correct transcript IDs for a gene
        """
        
        hgnc = "DYNLL1"
        lengths = _run(get_transcript_ids, hgnc)
        
        expected = ['ENST00000548342', 'ENST00000549989', 'ENST00000392509',
                    'ENST00000392508', 'ENST00000242577', 'ENST00000550845', 
                    'ENST00000550178', 'ENST00000548214', 'ENST00000552870',
                    'ENST00000549649']
        
        self.assertEqual(set(lengths), set(expected))
    
    def test_load_gene(self):
        """ check that we correctly load the suitable transcripts for a gene
        """
        
        # define a hgnc symbol to load transcript for, and de novo sites to
        # check against
        hgnc = "DYNLL1"
        sites = [120934226, 120936012]
        gene = _run(load_gene, hgnc)
        counts = minimise_transcripts(gene.transcripts, sites)
        
        # define the expected transcript, and make sure that it is in the list
        # of suitable transcripts. There can be multiple transcripts return if
        # more than one transcript of the maximal length includes all de novos.
        expected = self.set_transcript()
        self.assertIn(expected.get_name(), counts)
        
        # and make sure if none of the de novos fall in a suitable transcript,
        # then we get an empty dict.
        sites = [100, 200]
        gene = _run(load_gene, hgnc)
        transcripts = gene.transcripts
        counts = minimise_transcripts(transcripts, sites)
        self.assertEqual(counts, {})
    
    def test_count_de_novos_per_transcript(self):
        """ test that we count de novos in transcripts correctly
        """
        sites = [120934226, 120936012]
        
        expected = {'ENST00000549649': {'len': 126, 'n': 1},
            'ENST00000548214': {'len': 204, 'n': 1},
            'ENST00000242577': {'len': 270, 'n': 2},
            'ENST00000392508': {'len': 270, 'n': 2},
            'ENST00000392509': {'len': 270, 'n': 2},
            'ENST00000549989': {'len': 270, 'n': 2},
            'ENST00000550178': {'len': 204, 'n': 1},
            'ENST00000552870': {'len': 144, 'n': 1},
            'ENST00000550845': {'len': 204, 'n': 1},
            'ENST00000548342': {'len': 270, 'n': 2}}
        
        transcripts = [self._construct_transcript(x) for x in expected]
        
        counts = count_de_novos_per_transcript(transcripts, sites)
        self.assertEqual(counts, expected)
        
        # TODO: add test case for error from gene where no protein coding
        # TODO: transcript is available
        
    def test_minimise_transcripts(self):
        """ test that minimise_transcripts() works correctly
        """
        
        # run through a test case for a single gene
        sites = [120934226, 120936012]
        expected = {'ENST00000242577': {'len': 270, 'n': 2},
            'ENST00000392508': {'len': 270, 'n': 2},
            'ENST00000392509': {'len': 270, 'n': 2},
            'ENST00000549989': {'len': 270, 'n': 2},
            'ENST00000548342': {'len': 270, 'n': 2}}
        
        transcripts = [self._construct_transcript(x) for x in expected]
        counts = minimise_transcripts(transcripts, sites)
        
        self.assertEqual(counts, expected)
        
        # check that when we don't have any de novos, we return an empty list
        self.assertEqual(minimise_transcripts(transcripts, []), {})
        
        # check that when none of the de novos are in a transcript, we return
        # an empty list.
        self.assertEqual(minimise_transcripts(transcripts, [100]), {})
