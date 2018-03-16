""" class to test the EnsemblRequest class
"""

import json
import unittest
import time
import tempfile
import shutil

from denovonear.ensembl_requester import EnsemblRequest

class TestEnsemblRequestPy(unittest.TestCase):
    """ unit test the EnsemblRequest class
    """
    
    @classmethod
    def setUpClass(self):
        self.temp_dir = tempfile.mkdtemp()
        self.ensembl = EnsemblRequest(self.temp_dir, genome_build="grch37")
    
    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.temp_dir)
    
    def test_open_url(self):
        """ test that open_url() works correctly
        """
        
        headers = {"Content-Type": "application/json"}
        url = "http://rest.ensembl.org/overlap/id/ENSG00000172320?feature=gene"
        (response, status_code, headers) = self.ensembl.open_url(url, headers)
        
        response = json.loads(response)
        
        self.assertEqual(status_code, 200)
        self.assertEqual(response, [{
            "source": "ensembl_havana",
            "logic_name": "ensembl_havana_gene",
            "feature_type": "gene",
            "external_name": "OR5A1",
            "seq_region_name": "11",
            "strand": 1,
            "id": "ENSG00000172320",
            "gene_id": "ENSG00000172320",
            "version": 3,
            "assembly_name": "GRCh38",
            "description": "olfactory receptor family 5 subfamily A member 1 [Source:HGNC Symbol;Acc:HGNC:8319]",
            "end": 59451380,
            "biotype": "protein_coding",
            "start": 59436469}]
            )
    
    def test_get_genes_for_hgnc_id(self):
        """ test that get_genes_for_hgnc_id() works correctly
        """
        
        genes = self.ensembl.get_genes_for_hgnc_id("KMT2A")
        self.assertEqual(genes, ['ENSG00000118058', 'ENSG00000267910'])
    
    def test_get_previous_symbol(self):
        """ test that get_previous_symbol() works correctly
        """
        
        prev = self.ensembl.get_previous_symbol("KMT2A")
        self.assertEqual(prev, ["MLL"])
        
        # make a check for a gene with multiple documents, to check that we
        # don't raise an error
        prev = self.ensembl.get_previous_symbol("KRT16P1")
        self.assertEqual(prev, ["KRT14P"])
        
    def test_get_transcript_ids_for_ensembl_gene_ids(self):
        """ test that get_transcript_ids_for_ensembl_gene_ids() works correctly
        """
        
        hgnc = ["KMT2A", "MLL"]
        ensg = ['ENSG00000118058', 'ENSG00000267910']
        
        enst = self.ensembl.get_transcript_ids_for_ensembl_gene_ids(ensg, hgnc)
        
        self.assertEqual(set(enst), set(['ENST00000534358', 'ENST00000531904',
            'ENST00000389506', 'ENST00000354520', 'ENST00000532204',
            'ENST00000529852', 'ENST00000527869', 'ENST00000533790',
            'ENST00000392873']))
    
    def test_get_genomic_seq_for_transcript(self):
        """ check that get_genomic_seq_for_transcript() works correctly
        """
        
        seq = self.ensembl.get_genomic_seq_for_transcript("ENST00000302030", expand=0)
        
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
        
        seq = self.ensembl.get_cds_seq_for_transcript("ENST00000302030")
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
        
        seq = self.ensembl.get_protein_seq_for_transcript("ENST00000302030")
        self.assertEqual(seq, 'MSITKAWNSSSVTMFILLGFTDHPELQALLFVTFLGIYLTTLAWNLAL'
            'IFLIRGDTHLHTPMYFFLSNLSFIDICYSSAVAPNMLTDFFWEQKTISFVGCAAQFFFFVGMGLSE'
            'CLLLTAMAYDRYAAISSPLLYPTIMTQGLCTRMVVGAYVGGFLSSLIQASSIFRLHFCGPNIINHF'
            'FCDLPPVLALSCSDTFLSQVVNFLVVVTVGGTSFLQLLISYGYIVSAVLKIPSAEGRWKACNTCAS'
            'HLMVVTLLFGTALFVYLRPSSSYLLGRDKVVSVFYSLVIPMLNPLIYSLRNKEIKDALWKVLERKK'
            'VFS')
    
    def test_get_genomic_seq_for_region(self):
        """ test that get_genomic_seq_for_region() works correctly
        """
        
        # not that this test uses GRCh37 coordinates
        seq = self.ensembl.get_genomic_seq_for_region('11', 59210617, 59210637)
        self.assertEqual(seq, 'CTTGTCCTTGTGGTCCACGGG')
    
    def test_get_chrom_for_transcript(self):
        """ test that get_chrom_for_transcript() works correctly
        """
        
        chrom = self.ensembl.get_chrom_for_transcript("ENST00000534358", "KMT2A")
        self.assertEqual(chrom, "11")
    
    def test_get_exon_ranges_for_transcript(self):
        """ test that get_exon_ranges_for_transcript() works correctly
        """
        
        exons = self.ensembl.get_exon_ranges_for_transcript("ENST00000534358")
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
        
        cds = self.ensembl.get_cds_ranges_for_transcript("ENST00000534358")
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
    
    def test_rate_limit_ensembl_requests(self):
        """ test that rate_limit_ensembl_requests() works correctly
        """
        
        current_time = time.time()
        self.ensembl.prior_time = current_time
        
        self.ensembl.rate_limit_ensembl_requests()
        delta = self.ensembl.prior_time - current_time
        
        self.assertTrue(delta >= self.ensembl.rate_limit)
