""" load trincleotide mutation rates
"""

from __future__ import annotations

import warnings
from pkg_resources import resource_filename

__CYVCF2_ERROR = False
try:
    from cyvcf2 import VCF
except ImportError:
    __CYVCF2_ERROR = True

def load_mutation_rates(path=None):
    """ load sequence context-based mutation rates
    
    Args:
        path: path to table of sequence context-based mutation rates. If None,
            this defaults to per-trinucleotide rates provided by Kaitlin Samocha
            (Broad Institute).
    
    Returns:
        list of [initial, changed, rate] lists e.g. [['AGA', 'ATA', '5e-8']]
    """
    
    if path is None:
        path = resource_filename(__name__, "data/rates.txt")
    
    rates = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("from"): # ignore the header line
                continue
            
            line = [ x.encode('utf8') for x in line.strip().split() ]
            rates.append(line)
    
    return rates

class load_mutation_rates_in_region:
    ''' class to handle loading mutations rates within genome regions.

    This only permits loading data from files from Roulette e.g.
      http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/
    '''
    CHROMS = list(map(str, range(1, 23)))
    def __init__(self, paths: list[str]):
        ''' start with list of paths to all Roulette VCFs
        '''
        if __CYVCF2_ERROR:
            raise ValueError("this does not work due to problems installing cyvcf2")

        self.vcfs = {}
        for path in paths:
            tbx = VCF(path)
            # tag the tabixfile against all contigs (chromosomes) it contains
            for contig in self._check_chroms(tbx):
                if contig not in self.vcfs:
                    self.vcfs[contig] = []
                self.vcfs[contig].append(tbx)
    def _check_chroms(self, vcf):
        ''' find all the chromosomes actually present in a VCF
        '''
        contigs = []
        for chrom in self.CHROMS:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    next(vcf(f'{chrom}'))
                contigs.append(chrom)
            except StopIteration:
                pass
        return contigs
    def __call__(self, chrom: str, start: int, end: int, tag='MR') -> dict[int, dict[str, float]]:
        ''' get mutation rates within a genome region

        Args:
            chrom: chromosome
            start: start position of region to gather data for
            end: end position (inclusive) of region to gather data for
            tag: which rate type to extract data for. Options are:
                MR (Roulette) - the default.
                GR (GnomAD) - from Kaitlin Samocha's rates
                AR (Adjusted roulette) - unavailable for most sites
                MC (Carlsson et al 2018) - unavailable for many sites.
        '''
        # the Roulette VCFs lack 'chr' prefixes, change to match
        if chrom.startswith('chr'):
            chrom = chrom[3:]
        
        # run through all the VCFs that contain data for the given chromosome.
        # Storing the data in a dict should deduplicate sites, if necessary.
        data = {}
        for tbx in self.vcfs[chrom]:
            for row in tbx(f'{chrom}:{start}-{end+1}'):
                if row.POS not in data:
                    data[row.POS] = {}
                data[row.POS][row.ALT[0]] = row.INFO[tag]
        return data
