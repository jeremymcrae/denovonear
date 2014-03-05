""" function to translate DNA codons to single character amino acids
"""

from __future__ import division
from __future__ import print_function

import gzip
import os
import glob
import bisect


def get_filename_dict(folder):
    """ return a dict of filenames, indexed by their start chrom position
    
    Args:
        folder: folder for files for a single chromosome
    
    Returns:
        dict of filenames, indexed by their start chrom position
    """
    
    files = glob.glob(os.path.join(folder, "*phyloP46way.wigFix.gz"))
    
    folder_dict = {}
    for filename in files:
        name = os.path.basename(filename)
        start_pos = int(name.split(".")[1])
        folder_dict[start_pos] = name
    
    return folder_dict

def find_files(file_dict, start_pos, end_pos):
    """ finds the files that overlap the start and end positions
    
    Args:
        file_dict: dict of filenames, indexed by their start chrom position
        start_pos: chromosome position
        end_pos: chromosome position
    
    Returns:
        list of files that contain positions overlapping start_pos to end_pos
    """
    
    positions = sorted(file_dict)
    
    # find the positions of the first file, and the last file
    initial = bisect.bisect(positions, start_pos) - 1
    final = bisect.bisect(positions, end_pos) + 1
    
    files = []
    for pos in range(initial, final):
        files.append(file_dict[positions[pos]])
        if pos == len(positions) - 1: # don't go beyond the chrom end
            break
    
    return files

def load_conservation_file(folder, filename, start, end, scores):
    """ loads conservation scores from a wig file
    
    Args:
        folder: folder for files for a single chromosome
        filename: string of filename
        start: chromosome position
        end: chromosome position
        scores: dict of conservation scores, indexed by chromosome position
    
    Returns:
        scores, updated with additional posiotions from current file
    """
    
    path = os.path.join(folder, filename)
    
    with gzip.open(path, "rt") as f:
        pos = 0
        for line in f:
            pos += 1
            
            # occasionally the file jumps the position to a new site, not just
            # at the start of the file
            if line.startswith("fixed"):
                line = line.strip().split()
                pos = int(line[2].split("=")[1])
                continue
            
            if start < pos < end:  
                scores[pos] = float(line.strip())
            elif pos > end: # break out once we're past the gene
                break
    
    return scores

def load_conservation_scores(folder, chrom, start, end):
    """ load conservation scores for a chromosome region
    
    Args:
        folder: path to folder containing conservation score data
        chrom: chromosome number (1, 2 .. 22, X, Y)
        start: chromosome position at the start of the desired region
        end: chromosome position at the end of the desired region
    
    Returns:
        dict of conservation scores, indexed by chromosome position
    """
    
    # expand the boundaries a little, so we make sure everything is captured
    start -= 10
    end += 10
    
    chrom_folder = os.path.join(folder, "chr{0}".format(chrom))
    
    file_dict = get_filename_dict(chrom_folder)
    files = find_files(file_dict, start, end)
    
    scores = {}
    for filename in files:
        load_conservation_file(chrom_folder, filename, start, end, scores)
        
    return scores

