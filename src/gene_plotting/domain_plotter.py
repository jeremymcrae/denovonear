""" class to plot protein domains.
"""


class DomainPlotter(object):
    
    def plot_domains(self, domains, sequence, de_novos=None):
        """ plots protein domains
        """
        
        length = len(sequence) / 100
        
        # increment the y_offset position, so as to avoid overplotting between
        # gene, transcript and protein diagrams
        self.y_offset += self.box_height * 3
        
        # make sure the full protein is visible
        self.add_box(x_pos=0, width=100, color="white")
        
        for domain in domains:
            self.plot_single_domain(domain, length)
            
        for de_novo in de_novos:
            self.plot_single_denovo(de_novo, length)
    
    def plot_single_domain(self, domain, length):
        """ plots a single domain
        """
        
        x_pos = domain["start"] / length
        width = (domain["end"] - domain["start"]) / length
        x_center = x_pos + (width / 2)
        
        # add a box on the domain plot, as well as a text label centered 
        # below the box
        self.add_box(x_pos, width)
        self.add_text(x_center, domain["domain_type"], y_adjust=2, horizontalalignment="center")
    
    def plot_single_denovo(self, de_novo, length):
        """ plots a de novo on the domain plot
        """
        
        x_pos = de_novo / length
        width = 1
        
        self.add_de_novo(x_pos, width)


