
import gzip
from pathlib import Path
from typing import Dict, Tuple

def load_structure(path: Path) -> Dict[Tuple[str, int], dict[str, float]]:
    ''' a very basic PDB parser, to get positions of carbon atoms at the amino acids
    '''
    opener = gzip.open if str(path).endswith('gz') else open
    coords = {}
    with opener(path, 'rt') as handle:
        for line in handle:
            if not line.startswith('ATOM'):
                continue
            
            # only use the carbon atom position in each residue
            atom_name = line[12:16].strip(' ')
            if atom_name != 'C':
                continue
            
            residue = line[17:20]
            chain = line[21]
            residue_num = line[22:26]
            x_pos = float(line[30:38])
            y_pos = float(line[38:46])
            z_pos = float(line[46:54])
            
            key = (chain, residue_num)
            if key in coords:
                raise ValueError(f'seen this atom previously: {key}')
            coords[key] = {'residue': residue, 'x': x_pos, 'y': y_pos, 'z': z_pos}
    
    # TOOD: basic validation - do we have every residue included - check for gaps
    
    return coords
