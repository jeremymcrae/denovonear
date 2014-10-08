""" class to plot transcripts.
"""

from __future__ import division


class TranscriptPlotter(object):
    
    def plot_transcript(self, gene, de_novos=None):
        """ plots a transcript
        """
        
        # increment the y_offset position, so as to avoid overplotting between
        # gene, transcript and protein diagrams
        self.y_offset += self.box_height * 3
        
        strand = gene.strand
        length = (gene.get_coding_distance(gene.get_cds_start(), gene.get_cds_end()) + 1) / 100
        
        # give a label for the gene
        x_pos = 0
        self.add_text(x_pos, gene.get_name(), y_adjust=2)
        
        for start, end in gene.cds:
            x_pos = gene.get_coding_distance(gene.get_cds_start(), start) / length
            width = (gene.get_coding_distance(start, end) + 1) / length
            
            self.add_box(x_pos, width)
        


