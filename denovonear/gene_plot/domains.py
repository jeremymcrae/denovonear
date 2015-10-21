""" class to plot protein domains.
"""

from __future__ import division

class DomainPlot(object):
    """ class to plot protein domains.
    """
    
    def plot_domains(self, domains, sequence, de_novos=None):
        """ plots protein domains
        
        Args:
            domains: list of dictionaries for different domain annotations in a
                protein
            sequence: protein sequence as string
        """
        
        cds_coords = [ self.transcript.convert_chr_pos_to_cds_positions(x) for x in self.de_novos ]
        
        codons = [ self.transcript.get_codon_number_for_cds_position(x) for x in cds_coords ]
        intracodons = [ self.transcript.get_position_within_codon(x)/3 for x in cds_coords ]
        aa_coords = [ sum(x) for x in zip(codons, intracodons) ]
        
        length = len(sequence) / self.size
        
        # increment the y_offset position, so as to avoid overplotting between
        # gene, transcript and protein diagrams
        self.y_offset += self.box_height * 3
        
        # make sure the full protein is visible
        self.add_box(x_pos=0, width=self.size, fillcolor="white")
        
        for domain in domains:
            self.plot_single_domain(domain, length)
        
        coordinates = [ (x/length, 0.333/length) for x in aa_coords ]
        self.add_de_novos(coordinates)
        
        # and include the gene symbol and amino acid length on the domain plot
        text = "{} ({} aa)".format(self.hgnc_symbol, len(sequence))
        self.add_text(self.size, text, y_adjust=self.box_height*1.5, horizontalalignment="right")
    
    def plot_single_domain(self, domain, length):
        """ plots a single domain
        """
        
        x_pos = domain["start"] / length
        width = (domain["end"] - domain["start"]) / length
        x_center = x_pos + (width / 2)
        
        # add a box on the domain plot, as well as a text label centered
        # below the box
        self.add_box(x_pos, width, fillcolor="lightblue")
        self.add_text(x_center, domain["domain_type"], y_adjust=self.box_height*1.5, horizontalalignment="center")
       
