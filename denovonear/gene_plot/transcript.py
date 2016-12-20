""" class to plot transcripts.
"""

from __future__ import division


class TranscriptPlot(object):
    """ class to plot transcripts.
    """
    
    def plot_transcript(self):
        """ plots a transcript
        """
        
        # increment the y_offset position, so as to avoid overplotting between
        # gene, transcript and protein diagrams
        self.y_offset += self.box_height * 3
        
        cds_start = self.transcript.get_cds_start()
        cds_end = self.transcript.get_cds_end()
        length = (self.transcript.get_coding_distance(cds_start, cds_end) + 1) / self.size
        
        # give a label for the gene
        x_pos = 0
        self.add_text(x_pos, self.transcript.get_name(), y_adjust=self.box_height*1.5)
        
        for exon in self.transcript.get_cds():
            start, end = exon['start'], exon['end']
            if self.transcript.get_strand() == "+":
                x_pos = self.transcript.get_coding_distance(cds_start, start) / length
                width = (self.transcript.get_coding_distance(start, end) + 1) / length
            else:
                x_pos = self.transcript.get_coding_distance(cds_start, end) / length
                width = (self.transcript.get_coding_distance(start, end) + 1) / length
            
            self.add_box(x_pos, width, fillcolor="green")
        
        # and plot the de novo positions
        for x in self.de_novos:
            cds = self.transcript.chrom_pos_to_cds(self.de_novos[x]["start_pos"])
            x_pos = cds['pos']/length
            width = max(1/length, self.size/1000)
            
            self.de_novos[x]["coordinate"] = (x_pos, width)
            
        self.add_de_novos(self.de_novos)
        
        # and include the gene symbol and amino acid length on the domain plot
        text = "{} bp".format(self.transcript.get_coding_distance(cds_start, cds_end) + 1)
        self.add_text(self.size, text, y_adjust=self.box_height*1.5, horizontalalignment="right")
        
