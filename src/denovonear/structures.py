
import gzip
from gzip import GzipFile
import logging
from pathlib import Path
import tarfile
from typing import Dict, Iterable, List, Tuple

from gencodegenes.transcript import Transcript

def find_uniprot_matches(path: Path, uniprot_ids: List[str]) -> Iterable[GzipFile]:
    ''' search a tar file for pdb files matching uniprot IDs
    '''
    logging.info(f'checking for structures in {path} with {uniprot_ids}')
    with tarfile.open(path) as tarred:
        for x in tarred:
            if x.name.endswith('cif.gz'):
                continue
            for uniprot_id in uniprot_ids:
                if uniprot_id in x.name:
                    logging.info(f'found structure for {uniprot_id}')
                    yield gzip.open(tarred.extractfile(x.name), 'rt')

def load_structure(path: Path, uniprot_ids: List[str]) -> Dict[Tuple[str, int], Dict[str, float]]:
    ''' a very basic PDB parser, to get positions of carbon atoms at the amino acids
    '''
    coords = None
    for handle in find_uniprot_matches(path, uniprot_ids):
        if coords is not None:
            # we've hit a second file for the same gene, so skip out
            return None
        
        coords = {}
        for line in handle:
            if not line.startswith('ATOM'):
                continue
            
            # only use the carbon atom position in each residue
            atom_name = line[12:16].strip(' ')
            if atom_name != 'C':
                continue
            
            residue = line[17:20]
            chain = line[21]
            residue_num = int(line[22:26])
            x_pos = float(line[30:38])
            y_pos = float(line[38:46])
            z_pos = float(line[46:54])
            
            key = (chain, residue_num)
            if key in coords:
                raise ValueError(f'seen this atom previously: {key}')
            coords[key] = {'residue': residue, 'x': x_pos, 'y': y_pos, 'z': z_pos}
        
        # TOOD: basic validation - do we have every residue included - check for gaps
    
    return coords

def valid_structure(structure_coords: Dict[Tuple[str, int], Dict[str, float]],
                    transcript: Transcript) -> bool:
    ''' check the structure sufficiently matches the transcript
    '''
    
    # check there is only a single protein chain in the structure, fail if not
    n_chains = len(set(x[0] for x in structure_coords))
    if n_chains > 1:
        logging.warning(f'cannot get distances from structure with multiple chains')
        return False
    
    # make sure coords are linearly sorted without gaps
    residues = [x[1] for x in sorted(structure_coords)]
    if residues != list(range(1, len(structure_coords) + 1)):
        logging.warning(f'cannot get distances from structure - missing residues')
        return False
    
    # make sure the protein coords extend to the end of the transcript range
    last_residue = transcript.get_coding_distance(transcript.cds_end)['pos'] // 3
    if abs(last_residue - residues[-1]) > 2:
        logging.warning(f"transcript length doesn't match structure length")
        return False
    
    return True