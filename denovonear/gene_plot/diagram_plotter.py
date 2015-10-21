""" class to plot gene diagrams. This is a base class for plotting genes,
transcripts, and proteins (with domains).
"""

from __future__ import division
from __future__ import print_function

import webcolors
import cairocffi as cairo

from denovonear.gene_plot.domains import DomainPlot
from denovonear.gene_plot.transcript import TranscriptPlot
from denovonear.gene_plot.genomic import GenomicPlot

class DiagramPlotter(GenomicPlot, TranscriptPlot, DomainPlot):
    
    size = 1000
    y_offset = 0
    box_height = size/20
    
    def __init__(self, transcript, hgnc_symbol, de_novos, filename=None):
        """ start the class with a matplotlib plot
        
        Args:
            transcript: Transcript object for the gene
            hgnc_symbol: hgnc_symbol for the gene
            de_novos: list of genomic coordinates for de novo variants
        """
        
        self.transcript = transcript
        self.de_novos = de_novos
        self.hgnc_symbol = hgnc_symbol
        
        if filename is None:
            filename = "gene_plot.{0}".format(self.hgnc_symbol)
        
        self.surface = cairo.PDFSurface(filename + '.pdf',  self.size,  self.size/2)
        self.cr = cairo.Context(self.surface)
    
    def add_box(self, x_pos, width, y_adjust=0, height=None, **kwargs):
        """ adds a rectangle to the plot
        """
        
        if height == None:
            height = self.box_height
        
        if "strokecolor" in kwargs:
            strokecolor = webcolors.name_to_rgb(kwargs["strokecolor"], spec=u'css3')
            strokecolor = [ x/255 for x in strokecolor ]
        else:
            # default to a black stroke
            strokecolor = [0, 0, 0]
        
        if "fillcolor" in kwargs:
            fillcolor = webcolors.name_to_rgb(kwargs["fillcolor"], spec=u'css3')
            fillcolor = [ x/255 for x in fillcolor ]
        else:
            # default to a white fill
            fillcolor = [1, 1, 1]
        
        self.cr.set_line_width(self.size/200)
        self.cr.set_source_rgb(*strokecolor)
        self.cr.rectangle(x_pos, self.y_offset + y_adjust, width, height)
        self.cr.stroke_preserve()
        self.cr.set_source_rgb(*fillcolor)
        self.cr.fill()
    
    def add_text(self, x_pos, text, y_adjust=0, **kwargs):
        """ adds some text to the plot
        
        Args:
            x_pos: horizontal position at which to plot the text
            text: string of text to be plotted
            y_adjust: how far away to plot the text from the current y_offset
            kwargs: additional arguments, mainly to add arguments like
                horizontalalignment="center" to pyplot.text
        """
        
        y_pos = self.y_offset + y_adjust
        
        self.cr.select_font_face("Arial")
        self.cr.set_font_size(self.size/50)
        px = max(self.cr.device_to_user_distance(1, 1))
        fascent, fdescent, fheight, fxadvance, fyadvance = self.cr.font_extents()
        xbearing, ybearing, width, height, xadvance, yadvance = \
                self.cr.text_extents(text)
        
        if "horizontalalignment" in kwargs:
            if kwargs["horizontalalignment"] == "center":
                x_pos = x_pos - xbearing - width / 2
            if kwargs["horizontalalignment"] == "right":
                x_pos = x_pos - xbearing - width
        y_pos = y_pos - fdescent + fheight / 2
        
        self.cr.move_to(x_pos, y_pos)
        self.cr.set_source_rgb(0, 0, 0)
        self.cr.show_text(text)
    
    def add_de_novos(self, coordinates, color="red"):
        """ adds a line to show the position of a de novo
        """
        
        height = self.box_height / 2
        
        for i, (x_pos, width) in enumerate(coordinates):
            self.add_box(x_pos, width, height=height, y_adjust=-height, \
                fillcolor=color, strokecolor=color)
    
    def export_figure(self):
        """ exports the plot as a pdf
        """
        
        self.cr.save()
        
        self.surface.write_to_png("test.png")
        
