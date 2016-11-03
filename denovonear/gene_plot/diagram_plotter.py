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
from denovonear.gene_plot.consequences import Consequences

class DiagramPlotter(GenomicPlot, TranscriptPlot, DomainPlot, Consequences):
    
    # Define a set of colors that we can plot domains as, but don't allow any
    # colors close to white, as they are difficult to distingush from background
    global white
    white = ["aliceblue", "azure", "beige", "cornsilk", "floralwhite",
        "gainsboro", "ghostwhite", "honeydew", "ivory", "lavenderblush", "linen",
        "lightyellow", "mintcream", "oldlace", "seashell", "snow", "white",
        "whitesmoke"]
    colorset = webcolors.css3_names_to_hex.keys()
    colorset = [ x for x in colorset if x not in white ]
    
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
        
        if transcript.get_strand() == "-":
            for key in self.de_novos:
                self.de_novos[key]["start_pos"] -= 1
        
        if filename is None:
            filename = "gene_plot.{0}.pdf".format(self.hgnc_symbol)
        
        self.surface = cairo.PDFSurface(filename,  self.size,  self.size/1.5)
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
    
    def add_text(self, x_pos, text, y_adjust=0, horizontalalignment=None, rotate=None, fontsize=None, color=None):
        """ adds some text to the plot
        
        Args:
            x_pos: horizontal position at which to plot the text
            text: string of text to be plotted
            y_adjust: how far away to plot the text from the current y_offset
        """
        
        y_pos = self.y_offset + y_adjust
        
        if fontsize is None:
            fontsize = "medium"
        
        sizes = {"small": self.size/80, "medium": self.size/50, "large": self.size/30}
        
        self.cr.select_font_face("Arial")
        self.cr.set_font_size(sizes[fontsize])
        px = max(self.cr.device_to_user_distance(1, 1))
        fascent, fdescent, fheight, fxadvance, fyadvance = self.cr.font_extents()
        xbearing, ybearing, width, height, xadvance, yadvance = \
                self.cr.text_extents(text)
        
        if rotate is not None:
            x_rot_height = width * math.sin(math.radians(rotate))
            x_rot_width = math.sqrt(width ** 2 - x_rot_height ** 2)
            rot_x_delta = width - x_rot_width
            rot_y_delta = x_rot_height - height
        
        if horizontalalignment is not None:
            if horizontalalignment == "center":
                x_pos = x_pos - xbearing - width / 2
            elif horizontalalignment == "right":
                x_pos = x_pos - xbearing - width
            if rotate:
                x_pos = x_pos + rot_x_delta
                y_pos = y_pos - rot_y_delta
    
        y_pos = y_pos - fdescent + fheight / 2
        
        self.cr.move_to(x_pos, y_pos)
        
        if rotate is not None:
            self.cr.rotate(math.radians(rotate))
        
        if color is None:
            # default to a black fill
            color = [0, 0, 0]
        else:
            color = webcolors.name_to_rgb(color, spec=u'css3')
            color = [ x/255 for x in color ]
        
        self.cr.set_source_rgb(*color)
        self.cr.show_text(text)
        
        if rotate is not None:
            self.cr.rotate(-math.radians(rotate))
    
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
            
            self.add_box(*de_novos[key]["coordinate"], height=height, y_adjust=-height, \
                fillcolor=color, strokecolor=color, linewidth=0)
            
            (x, y) = self.rotate_indicator(de_novos, key, i, color, height)
                
            text = self.get_conseqence(de_novos[key]["start_pos"], \
                de_novos[key]["ref_allele"], de_novos[key]["alt_allele"])
            
            self.add_text(x, text, y, rotate=315, fontsize="small")
    
    def rotate_indicator(self, de_novos, key, i, color, height):
        """
        """
        
        (x_pos, width) = de_novos[key]["coordinate"]
        positions = self.get_overlap_positions(de_novos, key)
        if len(positions) == 1:
            return (x_pos, -height)
        
        radians = self.get_rotation(positions, i)
        dy = height * math.sin(radians)
        dx = math.sqrt(height**2 - dy**2)
        
        if math.degrees(radians) < 90:
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
        
        return (x_pos + dx, -height-dy)
    
    def get_overlap_positions(self, de_novos, key):
        """ identify the variants which overlap the current site
        
        Args:
            de_novos: list of de novos within gene
            key: current de novo to check
        
        Returns:
            list positions for overlapping variants
        """
        
        (x_pos, width) = de_novos[key]["coordinate"]
        
        # check how many other sites the de novo overlaps
        overlaps = [ x_pos <= de_novos[z]["coordinate"][0] + de_novos[z]["coordinate"][1] \
            and x_pos + width >= de_novos[z]["coordinate"][0] for z in sorted(de_novos) ]
        return [ x for x, value in enumerate(overlaps) if value ]
    
    def get_rotation(self, positions, i):
        """ get the rotation angle for line segments for overlapping de novos
        
        Args:
            positions: list of positions that overlap the current de novo
            i: position in list for current de novo
        
        Returns:
            angle to plot segment in radians
        """
        
        # figure out the angle for the indicator line
        increment = 180/(len(positions) + 1)
        angle = increment + positions.index(i) * increment
        return math.radians(angle)
    
    def export_figure(self):
        """ exports the plot as a pdf
        """
        
        self.cr.save()
        
