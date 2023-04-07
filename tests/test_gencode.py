
import asyncio
from pathlib import Path
import unittest
import tempfile

from denovonear.rate_limiter import RateLimiter
from denovonear.load_gene import construct_gene_object
from gencodegenes.gencode import Gencode, _parse_gtfline, _open_gencode

async def call(func, *args, **kwargs):
    ''' call ensembl rest API function
    '''
    async with RateLimiter(15) as ensembl:
        return await func(ensembl, *args, **kwargs)

def _run(func, *args, **kwargs):
    return asyncio.run(call(func, *args, **kwargs))

def write_gtf(path, lines):
    with open(path, 'wt') as output:
        for line in lines:
            output.write(line)

def make_fasta(path, chroms):
    with open(path, 'wt') as output:
        for chrom in chroms:
            output.write(f'>{chrom}\n')
            lines = ['A' * 50 + '\n'] * 50
            output.writelines(lines)

class TestGencode(unittest.TestCase):
    
    def setUp(self):
        ''' set path to folder with test data
        '''
        self.folder = Path(__file__).parent.parent /  "data"
        self.gtf_path = self.folder / 'example.grch38.gtf'
        self.fasta_path = self.folder / 'example.grch38.fa'
        temp_gtf = tempfile.NamedTemporaryFile(delete=False)
        temp_fasta = tempfile.NamedTemporaryFile(delete=False)
        self.temp_gtf_path = temp_gtf.name
        self.temp_fasta_path = temp_fasta.name
        temp_gtf.close()
        temp_fasta.close()
        self.maxDiff = None
    
    def tearDown(self):
        try:
            Path(self.temp_gtf_path).unlink()
            Path(self.temp_fasta_path).unlink()
        except:
            pass
        
    def test_gencode_matches_ensembl(self):
        ''' test thaqt gencode data matches ensembl
        '''
        gencode = Gencode(self.gtf_path, self.fasta_path)
        for gencode_tx in gencode['OR4F5'].transcripts:
            # get the transcript ID (but trim the version number)
            tx_id = gencode_tx.name.split('.')[0]
            ensembl_tx = _run(construct_gene_object, tx_id)
            self.assertEqual(gencode_tx.exons, ensembl_tx.exons)
            self.assertEqual(gencode_tx.cds, ensembl_tx.cds)
            self.assertEqual(gencode_tx.cds_sequence, ensembl_tx.cds_sequence)
