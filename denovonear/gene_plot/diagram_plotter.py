""" class to plot gene diagrams. This is a base class for plotting genes,
transcripts, and proteins (with domains).
"""

from __future__ import division
from __future__ import print_function

import math

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
        
        self.surface = cairo.PDFSurface(filename + '.pdf',  self.size,  self.size/1.5)
        self.cr = cairo.Context(self.surface)
    
    def add_box(self, x_pos, width, y_adjust=0, height=None, strokecolor=None, fillcolor=None, linewidth=None):
        """ adds a rectangle to the plot
        """
        
        if height is None:
            height = self.box_height
        
        if strokecolor is None:
            # default to a black stroke
            strokecolor = [0, 0, 0]
        else:
            strokecolor = webcolors.name_to_rgb(strokecolor, spec=u'css3')
            strokecolor = [ x/255 for x in strokecolor ]
        
        if fillcolor is None:
            # default to a white fill
            fillcolor = [1, 1, 1]
        else:
            fillcolor = webcolors.name_to_rgb(fillcolor, spec=u'css3')
            fillcolor = [ x/255 for x in fillcolor ]
        
        if linewidth is None:
            linewidth = self.size/200
        
        self.cr.set_line_width(linewidth)
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
    
    def add_de_novos(self, de_novos):
        """ adds a line to show the position of a de novo
        
        Args:
            de_novos: dictionary of variants, indexed by (position, person_id)
                tuples. Each variant entry contains a "coordinate" entry, which
                is a tuple of (x_position, width).
        """
        
        height = self.box_height / 2
        
        i = -1
        for key in sorted(de_novos):
            i += 1
            color = "gray"
            if de_novos[key]["source"] == "internal":
                color = "red"
            
            coords = de_novos[key]["coordinate"]
            x_pos = coords[0]
            width = coords[1]
            
            self.add_box(x_pos, width, height=height, y_adjust=-height, \
                fillcolor=color, strokecolor=color, linewidth=0)
            
            # check how many other sites the de novo overlaps
            overlaps = [ x_pos <= de_novos[z]["coordinate"][0] + de_novos[z]["coordinate"][1] \
                and x_pos + width >= de_novos[z]["coordinate"][0] for z in sorted(de_novos) ]
            positions = [ x for x, value in enumerate(overlaps) if value ]
            n_overlaps = sum(overlaps)
            
            if n_overlaps == 1:
                continue
            
            # figure out the angle for the indicator line
            increment = 180/(n_overlaps + 1)
            angle = increment + positions.index(i) * increment
            radians = math.radians(angle)
            
            dy = height * math.sin(radians)
            dx = math.sqrt(height**2 - dy**2)
            
            if angle < 90:
                dx = -dx
            
            # get the rgb values for the required stroke color
            strokecolor = webcolors.name_to_rgb(color, spec=u'css3')
            strokecolor = [ x/255 for x in strokecolor ]

            # plot a line arcing out from the overlap point
            self.cr.new_path()
            self.cr.move_to(x_pos, self.y_offset - height)
            self.cr.rel_line_to(dx, -dy)
            self.cr.set_line_width(self.size/1000)
            self.cr.set_source_rgb(*strokecolor)
            self.cr.stroke()
    
    def export_figure(self):
        """ exports the plot as a pdf
        """
        
        self.cr.save()
        
        self.surface.write_to_png("test.png")
        
