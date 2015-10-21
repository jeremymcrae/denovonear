""" class to plot transcripts.
"""

from __future__ import division


class TranscriptPlot(object):
    """ class to plot transcripts.
    """
    
    def plot_transcript(self):
        """ plots a transcript
        """
        
        cds_coords = [ self.transcript.convert_chr_pos_to_cds_positions(x) for x in self.de_novos ]
        
        # increment the y_offset position, so as to avoid overplotting between
        # gene, transcript and protein diagrams
        self.y_offset += self.box_height * 3
        
        cds_start = self.transcript.get_cds_start()
        cds_end = self.transcript.get_cds_end()
        length = (self.transcript.get_coding_distance(cds_start, cds_end) + 1) / self.size
        
        # give a label for the gene
        x_pos = 0
        self.add_text(x_pos, self.transcript.get_name(), y_adjust=self.box_height*1.5)
        
        for start, end in self.transcript.cds:
            if self.transcript.strand == "+":
                x_pos = self.transcript.get_coding_distance(cds_start, start) / length
                width = (self.transcript.get_coding_distance(start, end) + 1) / length
            else:
                x_pos = self.transcript.get_coding_distance(cds_start, end) / length
                width = (self.transcript.get_coding_distance(start, end) + 1) / length
            
            self.add_box(x_pos, width, fillcolor="green")
        
        # and plot the de novo positions
        coordinates = [ (x/length, 1/length) for x in cds_coords ]
        self.add_de_novos(coordinates)
        
