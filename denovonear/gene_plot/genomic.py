""" class to plot genes.
"""

from __future__ import division


class GenomicPlot(object):
    """ class to plot exons of genes using genomic coordinates.
    """
    
    def plot_gene(self):
        """ plots the exons of a gene
        """
        
        # increment the y_offset position, so as to avoid overplotting between
        # gene, transcript and protein diagrams
        self.y_offset += self.box_height * 3
        
        min_pos = self.transcript.get_start()
        max_pos = self.transcript.get_end()
        strand = self.transcript.get_strand()
        length = (max_pos - min_pos) / self.size
        
        # plot the base strand
        x_pos = (self.transcript.get_start() - min_pos) / length
        width = (self.transcript.get_end() - self.transcript.get_start()) / length
        height = self.box_height/8
        y_adjust = self.box_height/2 - height/2
        self.add_box(x_pos, width, height=height, y_adjust=y_adjust, fillcolor="black")
        
        # give a label for the gene
        self.add_text(x_pos, self.transcript.get_name(), y_adjust=self.box_height*1.5)
        
        for exon in self.transcript.get_exons():
            self.plot_single_exon(exon['start'], exon['end'], length, min_pos)
        
        for x in self.de_novos:
            x_pos = (self.de_novos[x]["start_pos"] - min_pos)/length
            width = max(1/length, self.size/1000)
            self.de_novos[x]["coordinate"] = (x_pos, width)
        
        self.add_de_novos(self.de_novos)
    
    def mixed_coords(self, start, end, length, min_pos):
        """
        """
        
        cds_min = min(self.transcript.get_cds_start(), self.transcript.get_cds_end())
        cds_max = max(self.transcript.get_cds_start(), self.transcript.get_cds_end())
        strand = self.transcript.get_strand()
        
        x_pos_1 = (start - min_pos) / length
        width_1 = (end - start) / length
        
        if start <= cds_min <= end:
            midpoint = cds_min
            color_1 = "white"
            color_2 = "green"
        else:
            midpoint = cds_max
            color_1 = "green"
            color_2 = "white"
        
        x_pos_2 = (midpoint - min_pos) / length
        width_2 = (end - midpoint) / length
        
        return (x_pos_1, width_1, color_1), (x_pos_2, width_2, color_2)
        
    def plot_single_exon(self, start, end, length, min_pos):
        """ adds a rectangle to the plot
        
        NOTE: this currently doesn't allow for single exon genes where the
        single exon has 5' and 3' untranslated regions as well as a coding
        region within the exon.
        """
        
        cds_min = min(self.transcript.get_cds_start(), self.transcript.get_cds_end())
        cds_max = max(self.transcript.get_cds_start(), self.transcript.get_cds_end())
        
        x_pos = (start - min_pos) / length
        width = (end - start) / length
        
        if self.transcript.in_coding_region(start) and self.transcript.in_coding_region(end):
            self.add_box(x_pos, width, fillcolor="green")
        elif not self.transcript.in_coding_region(start) and not self.transcript.in_coding_region(end):
            self.add_box(x_pos, width, fillcolor="white")
        else:
            # exons with coding and untranslated regions, either utr first, or
            # utr second
            (x_pos_1, width_1, color_1), (x_pos_2, width_2, color_2) = self.mixed_coords(start, end, length, min_pos)
            
            self.add_box(x_pos_1, width_1, fillcolor=color_1)
            self.add_box(x_pos_2, width_2, fillcolor=color_2)
