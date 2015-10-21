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
        
        length = (gene.get_coding_distance(gene.get_cds_start(), gene.get_cds_end()) + 1) / self.size
        
        # give a label for the gene
        x_pos = 0
        self.add_text(x_pos, gene.get_name(), y_adjust=self.box_height*1.5)
        
        for start, end in gene.cds:
            if gene.strand == "+":
                x_pos = gene.get_coding_distance(gene.get_cds_start(), start) / length
                width = (gene.get_coding_distance(start, end) + 1) / length
            else:
                x_pos = gene.get_coding_distance(gene.get_cds_start(), end) / length
                width = (gene.get_coding_distance(start, end) + 1) / length
            
            self.add_box(x_pos, width, fillcolor="green")
        
        # and plot the de novo positions
        for de_novo in de_novos:
            x_pos = de_novo / length
            width = 1/length
            self.add_de_novo(x_pos, width)
