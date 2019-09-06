""" class to test the EnsemblRequest class
"""

import json
import os
from pathlib import Path
import asyncio
import unittest
import time
import tempfile
import shutil

from denovonear.rate_limiter import RateLimiter
from denovonear.ensembl_requester import (get_genes_for_hgnc_id,
    get_transcript_ids_for_ensembl_gene_id, get_previous_symbol,
    get_genomic_seq_for_transcript, get_cds_seq_for_transcript,
    get_protein_seq_for_transcript, get_exon_ranges_for_transcript,
    get_cds_ranges_for_transcript)

async def call(func, *args, **kwargs):
    ''' call ensembl rest API function
    '''
    async with RateLimiter(15) as ensembl:
        return await func(ensembl, *args, **kwargs)

def _run(func, *args, **kwargs):
    return asyncio.get_event_loop().run_until_complete(call(func, *args, **kwargs))

class TestEnsemblRequestPy(unittest.TestCase):
    """ unit test the EnsemblRequest class
    """
    
    def test_get_genes_for_hgnc_id(self):
        """ test that get_genes_for_hgnc_id() works correctly
        """
        
        genes = _run(get_genes_for_hgnc_id, "KMT2A", build='grch37')
        self.assertEqual(genes, ['ENSG00000118058', 'ENSG00000267910'])
    
    def test_get_previous_symbol(self):
        """ test that get_previous_symbol() works correctly
        """
        
        prev = _run(get_previous_symbol, "KMT2A")
        self.assertEqual(prev, ["MLL"])
        
        # make a check for a gene with multiple documents, to check that we
        # don't raise an error
        prev = _run(get_previous_symbol, "KRT16P1")
        self.assertEqual(prev, ["KRT14P"])
    
    def test_get_transcript_ids_for_ensembl_gene_id(self):
        """ test that get_transcript_ids_for_ensembl_gene_id() works correctly
        """
        
        hgnc = ["KMT2A", "MLL"]
        
        enst = _run(get_transcript_ids_for_ensembl_gene_id, 'ENSG00000118058', hgnc)
        enst += _run(get_transcript_ids_for_ensembl_gene_id, 'ENSG00000267910', hgnc)
        
        self.assertEqual(set(enst), set(['ENST00000534358', 'ENST00000531904',
            'ENST00000389506', 'ENST00000354520', 'ENST00000532204',
            'ENST00000529852', 'ENST00000527869', 'ENST00000533790',
            'ENST00000392873']))
    
    def test_get_genomic_seq_for_transcript(self):
        """ check that get_genomic_seq_for_transcript() works correctly
        """
        
        seq = _run(get_genomic_seq_for_transcript, "ENST00000302030", expand=0)
        
        self.assertEqual(seq, ('11', 59210617, 59211667, '+', 'CTTGTCCTTGTGGTCC'
            'ACGGGAAGCATGTCCATAACCAAAGCCTGGAACAGCTCATCAGTGACCATGTTCATCCTCCTGGGA'
            'TTCACAGACCATCCAGAACTCCAGGCCCTCCTCTTTGTGACCTTCCTGGGCATCTATCTTACCACC'
            'CTGGCCTGGAACCTGGCCCTCATTTTTCTGATCAGAGGTGACACCCATCTGCACACACCCATGTAC'
            'TTCTTCCTAAGCAACTTATCTTTCATTGACATCTGCTACTCTTCTGCTGTGGCTCCCAATATGCTC'
            'ACTGACTTCTTCTGGGAGCAGAAGACCATATCATTTGTGGGCTGTGCTGCTCAGTTTTTTTTCTTT'
            'GTCGGCATGGGTCTGTCTGAGTGCCTCCTCCTGACTGCTATGGCATACGACCGATATGCAGCCATC'
            'TCCAGCCCCCTTCTCTACCCCACTATCATGACCCAGGGCCTCTGTACACGCATGGTGGTTGGGGCA'
            'TATGTTGGTGGCTTCCTGAGCTCCCTGATCCAGGCCAGCTCCATATTTAGGCTTCACTTTTGCGGA'
            'CCCAACATCATCAACCACTTCTTCTGCGACCTCCCACCAGTCCTGGCTCTGTCTTGCTCTGACACC'
            'TTCCTCAGTCAAGTGGTGAATTTCCTCGTGGTGGTCACTGTCGGAGGAACATCGTTCCTCCAACTC'
            'CTTATCTCCTATGGTTACATAGTGTCTGCGGTCCTGAAGATCCCTTCAGCAGAGGGCCGATGGAAA'
            'GCCTGCAACACGTGTGCCTCGCATCTGATGGTGGTGACTCTGCTGTTTGGGACAGCCCTTTTCGTG'
            'TACTTGCGACCCAGCTCCAGCTACTTGCTAGGCAGGGACAAGGTGGTGTCTGTTTTCTATTCATTG'
            'GTGATCCCCATGCTGAACCCTCTCATTTACAGTTTGAGGAACAAAGAGATCAAGGATGCCCTGTGG'
            'AAGGTGTTGGAAAGGAAGAAAGTGTTTTCTTAGGTCATGCGTAGAAACTTATTTATCCAAACTGCT'
            'GGAGAATTAAACAATCCAAGCCTTCACCTCCACCTCTGCCTCAGG'))
    
    def test_get_cds_seq_for_transcript(self):
        """ check that get_cds_seq_for_transcript() works correctly
        """
        
        seq = _run(get_cds_seq_for_transcript, "ENST00000302030")
        self.assertEqual(seq, 'ATGTCCATAACCAAAGCCTGGAACAGCTCATCAGTGACCATGTTCATC'
            'CTCCTGGGATTCACAGACCATCCAGAACTCCAGGCCCTCCTCTTTGTGACCTTCCTGGGCATCTAT'
            'CTTACCACCCTGGCCTGGAACCTGGCCCTCATTTTTCTGATCAGAGGTGACACCCATCTGCACACA'
            'CCCATGTACTTCTTCCTAAGCAACTTATCTTTCATTGACATCTGCTACTCTTCTGCTGTGGCTCCC'
            'AATATGCTCACTGACTTCTTCTGGGAGCAGAAGACCATATCATTTGTGGGCTGTGCTGCTCAGTTT'
            'TTTTTCTTTGTCGGCATGGGTCTGTCTGAGTGCCTCCTCCTGACTGCTATGGCATACGACCGATAT'
            'GCAGCCATCTCCAGCCCCCTTCTCTACCCCACTATCATGACCCAGGGCCTCTGTACACGCATGGTG'
            'GTTGGGGCATATGTTGGTGGCTTCCTGAGCTCCCTGATCCAGGCCAGCTCCATATTTAGGCTTCAC'
            'TTTTGCGGACCCAACATCATCAACCACTTCTTCTGCGACCTCCCACCAGTCCTGGCTCTGTCTTGC'
            'TCTGACACCTTCCTCAGTCAAGTGGTGAATTTCCTCGTGGTGGTCACTGTCGGAGGAACATCGTTC'
            'CTCCAACTCCTTATCTCCTATGGTTACATAGTGTCTGCGGTCCTGAAGATCCCTTCAGCAGAGGGC'
            'CGATGGAAAGCCTGCAACACGTGTGCCTCGCATCTGATGGTGGTGACTCTGCTGTTTGGGACAGCC'
            'CTTTTCGTGTACTTGCGACCCAGCTCCAGCTACTTGCTAGGCAGGGACAAGGTGGTGTCTGTTTTC'
            'TATTCATTGGTGATCCCCATGCTGAACCCTCTCATTTACAGTTTGAGGAACAAAGAGATCAAGGAT'
            'GCCCTGTGGAAGGTGTTGGAAAGGAAGAAAGTGTTTTCTTAG')
    
    def test_get_protein_seq_for_transcript(self):
        """ test that get_protein_seq_for_transcript() works correctly
        """
        
        seq = _run(get_protein_seq_for_transcript, "ENST00000302030")
        self.assertEqual(seq, 'MSITKAWNSSSVTMFILLGFTDHPELQALLFVTFLGIYLTTLAWNLAL'
            'IFLIRGDTHLHTPMYFFLSNLSFIDICYSSAVAPNMLTDFFWEQKTISFVGCAAQFFFFVGMGLSE'
            'CLLLTAMAYDRYAAISSPLLYPTIMTQGLCTRMVVGAYVGGFLSSLIQASSIFRLHFCGPNIINHF'
            'FCDLPPVLALSCSDTFLSQVVNFLVVVTVGGTSFLQLLISYGYIVSAVLKIPSAEGRWKACNTCAS'
            'HLMVVTLLFGTALFVYLRPSSSYLLGRDKVVSVFYSLVIPMLNPLIYSLRNKEIKDALWKVLERKK'
            'VFS')
    
    def test_get_exon_ranges_for_transcript(self):
        """ test that get_exon_ranges_for_transcript() works correctly
        """
        
        exons = _run(get_exon_ranges_for_transcript, "ENST00000534358")
        self.assertEqual(exons, [(118307205, 118307659), (118339490, 118339559),
            (118342377, 118345030), (118347520, 118347697), (118348682, 118348916),
            (118350889, 118350953), (118352430, 118352807), (118353137, 118353210),
            (118354898, 118355029), (118355577, 118355690), (118359329, 118359475),
            (118360507, 118360602), (118360844, 118360964), (118361911, 118362033),
            (118362459, 118362643), (118363772, 118363945), (118365003, 118365113),
            (118365409, 118365482), (118366415, 118366608), (118366976, 118367082),
            (118368651, 118368788), (118369085, 118369243), (118370018, 118370135),
            (118370550, 118370628), (118371702, 118371862), (118372387, 118372572),
            (118373113, 118377361), (118378244, 118378324), (118379851, 118379915),
            (118380663, 118380833), (118382666, 118382740), (118390333, 118390507),
            (118390672, 118390779), (118391517, 118391600), (118392003, 118392132),
            (118392612, 118397539)])
    
    def test_get_cds_ranges_for_transcript(self):
        """ tets that get_cds_ranges_for_transcript() works correctly
        """
        
        cds = _run(get_cds_ranges_for_transcript, "ENST00000534358")
        self.assertEqual(cds, [(118307228, 118307659), (118339490, 118339559),
            (118342377, 118345030), (118347520, 118347697), (118348682, 118348916),
            (118350889, 118350953), (118352430, 118352807), (118353137, 118353210),
            (118354898, 118355029), (118355577, 118355690), (118359329, 118359475),
            (118360507, 118360602), (118360844, 118360964), (118361911, 118362033),
            (118362459, 118362643), (118363772, 118363945), (118365003, 118365113),
            (118365409, 118365482), (118366415, 118366608), (118366976, 118367082),
            (118368651, 118368788), (118369085, 118369243), (118370018, 118370135),
            (118370550, 118370628), (118371702, 118371862), (118372387, 118372572),
            (118373113, 118377361), (118378244, 118378324), (118379851, 118379915),
            (118380663, 118380833), (118382666, 118382740), (118390333, 118390507),
            (118390672, 118390779), (118391517, 118391600), (118392003, 118392132),
            (118392612, 118392887)])
