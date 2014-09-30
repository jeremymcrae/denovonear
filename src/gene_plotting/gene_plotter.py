""" class to plot genes.
"""

from __future__ import division


class GenePlotter(object):
    
    def plot_gene(self, gene, min_pos, max_pos, de_novos=None):
        """ adds a rectangle to the plot
        """
        
        # increment the y_offset position, so as to avoid overplotting between
        # gene, transcript and protein diagrams
        self.y_offset += self.box_height * 3
        
        cds_min = min(gene.get_cds_start(), gene.get_cds_end())
        cds_max = max(gene.get_cds_start(), gene.get_cds_end())
        
        strand = gene.strand
        length = (max_pos - min_pos) / 100
        
        # plot the base strand
        x_pos = (gene.get_start() - min_pos) / length
        width = (gene.get_end() - gene.get_start()) / length
        height = self.box_height/8
        y_adjust = self.box_height/2 - height/2
        self.add_box(x_pos, width, height=height, y_adjust=y_adjust, alpha=1, color="black")
        
        # give a label for the gene
        self.add_text(x_pos, gene.get_name(), y_adjust=2)
        
        for (start, end) in gene.exons:
            self.plot_single_exon(start, end, length, gene, cds_min, cds_max)
        
    def plot_single_exon(self, start, end, length, gene, cds_min, cds_max):
        x_pos = (start - min_pos) / length
        width = (end - start) / length
        
        if gene.in_coding_region(start) and gene.in_coding_region(end):
            self.add_box(x_pos, width)
        elif not gene.in_coding_region(start) and not gene.in_coding_region(end):
            self.add_box(x_pos, width, color="white")
        elif not gene.in_coding_region(start) and gene.in_coding_region(end):
            width = (cds_min - start) / length
            self.add_box(x_pos, width, color="white")
            x_pos = cds_min - min_pos / length
            width = (end - cds_min) / length
            self.add_box(x_pos, width)
        elif gene.in_coding_region(start) and not gene.in_coding_region(end):
            width = (cds_max - start) / length
            self.add_box(x_pos, width)
            x_pos = cds_max - min_pos / length
            width = (end - cds_max) / length
            self.add_box(x_pos, width, color="white")
