""" function to translate DNA codons to single character amino acids
"""

from __future__ import division
from __future__ import print_function

import os
import gzip

CONSERVATION_FOLDER = "/lustre/scratch109/sanger/jm33/phyloP46way/vertebrate/"

def export_gzip(folder, chrom, new_file):
    """ write the previous lines to a gzip file
    """
    
    new_pos = new_file[0].split()[2].split("=")[1]
    new_file = "".join(new_file)
    new_path = "chr{0}.{1}.phyloP46way.wigFix.gz".format(chrom, new_pos)
    new_path = os.path.join(folder, new_path)
    with gzip.open(new_path, 'wt') as f:
        f.write(new_file)

def split_phyloP_into_intervals(folder, chrom):
    """ splits phyloP wig files into files containing different intervals, so
    we can quickly load a single file containing the positions we want, rather 
    than having to load a full chromosome file for a single rehion.
    """
    
    filename = "chr{0}.phyloP46way.wigFix.gz".format(chrom)
    path = os.path.join(folder, filename)
    
    # figure out the folder name for the chromosome specific files
    new_folder = os.path.join(folder, "chr{0}".format(chrom))
    if not os.path.exists(new_folder):
        os.mkdir(new_folder)
    
    # run through the old chromosome file
    with gzip.open(path, "rt") as f:
        start_pos = None
        pos = 0
        new_file = []
        for line in f:
            pos += 1
            
            if start_pos is None:
                start_pos = pos
                new_file.append(line)
            elif line.startswith("fixed") and abs(start_pos - pos) > 500000:
                tmp = line.strip().split()
                pos = int(tmp[2].split("=")[1])
                print(new_file[0].strip())
                new_pos = new_file[0].split()[2].split("=")[1]
                
                # write the previous lines to a gzip file
                export_gzip(new_folder, chrom, new_file)
                
                start_pos = pos
                new_file = [line]
            else:
                new_file.append(line)

def main():
    """
    """
    
    chroms = list(range(1, 22))
    # chroms = list(range(9, 22))
    chroms.append("X")
    chroms.append("Y")
    
    # split_phyloP_into_intervals(CONSERVATION_FOLDER, 22)
    
    for chrom in chroms:
        split_phyloP_into_intervals(CONSERVATION_FOLDER, chrom)

if __name__ == '__main__':
    main()





