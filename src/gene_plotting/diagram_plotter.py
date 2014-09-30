""" class to plot gene diagrams. This is a base class for plotting genes,
transcripts, and proteins (with domains).
"""

from __future__ import division
from __future__ import print_function


import matplotlib
matplotlib.use("pdf")
from matplotlib import pyplot

from src.gene_plotting.domain_plotter import DomainPlotter
from src.gene_plotting.transcript_plotter import TranscriptPlotter
from src.gene_plotting.gene_plotter import GenePlotter

class DiagramPlotter(GenePlotter, TranscriptPlotter, DomainPlotter):
    
    y_offset = 0
    box_height = 5
    
    def __init__(self):
        """ start the class with a matplotlib plot
        """
        
        figure = pyplot.figure()
        self.plot = figure.add_subplot(1, 1, 1)
        self.de_novo_positions = []
    
    def add_box(self, x_pos, width, color="green", alpha=0.5, height=None, y_adjust=0):
        """ adds a rectangle to the plot
        """
        
        if height == None:
            height = self.box_height
        
        xy = (x_pos, self.y_offset + y_adjust)
        rect = matplotlib.patches.Rectangle(xy, width, height=height)
        rect.set_facecolor(color)
        rect.set_alpha(alpha)
        
        self.plot.add_patch(rect)
    
    def add_text(self, x_pos, text, y_adjust=0, **kwargs):
        """ adds some text to the plot
        
        Args:
            x_pos: horizontal position at which to plot the text
            text: string of text to be plotted
            y_adjust: how far away to plot the text from the current y_offset
            kwargs: additional arguments, mainly to add arguments like 
                horizontalalignment="center" to pyplot.text
        """
        
        y_pos = self.y_offset - y_adjust
        xy = (x_pos, y_pos)
        
        pyplot.text(x_pos, y_pos, text, size="x-large", family="sans-serif", **kwargs)
    
    def add_de_novo(self, x_pos, width):
        """ adds a line to show the position of a de novo
        """
        
        self.de_novo_positions.append(x_pos)
        
        height = self.box_height / 2
        y_adjust = -self.box_height
        
        self.add_box(x_pos, width, color="red", height=height, y_adjust=y_adjust)
    
    def export_figure(self, path="test.pdf"):
        """ exports the plot as a pdf
        """
        
        # scale the axes, then remove them from view
        pyplot.axis("off")
        self.plot.autoscale()
        # self.plot.autoscale_view()
        
        # and save the plot to a file
        pyplot.savefig(path, bbox_inches="tight")
        

