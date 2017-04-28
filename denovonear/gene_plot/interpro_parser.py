""" functions to parse the output from requests to the InterPro REST API
"""

# from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys

class InterProParser(object):
    """ class to parse results from domain prediction runs from the EBI REST API
    """
    
    def find_domain_match(self, entry, domains):
        """ find matches between the current domain to the other domain results
        
        If the domain overlaps any other domain, we check if the description
        terms match, if they do, it is simply a matter of different tools
        predicting the same domain, in which case we avoid duplicating the
        domain.
        
        Args:
            entry: dictionary entry for domain prediction
        
        Returns:
            Return -1 when there are no matches, otherwise return the list
            position that the match is found at
        """
        
        match = -1
        for index in range(len(domains)):
            result = domains[index]
            if entry["end"] >= result["start"] and result["end"] >= entry["start"]:
                if len(entry["domain_type"] & result["domain_type"]) > 0:
                    match = index
                    break
        
        return match
    
    def standardise_names(self, domains):
        """ picks the shortest name from a list of alternates for each domain
        """
        
        # convert the list of alternative domain names into a single identifier
        for pos in range(len(domains)):
            names = list(domains[pos]["domain_type"])
            
            names = [ x.replace(" domain", "") for x in names ]
            name_lengths = [ len(x) for x in names ]
            
            # find the shortest name (note, this won't necessarily be the best
            # name, but it's a start).
            name_pos = name_lengths.index(min(name_lengths))
            
            domains[pos]["domain_type"] = names[name_pos]
        
        return domains
    
    def build_domain_entry(self, line):
        """ clean up a interpro results line.
        """
        
        entry = {}
            
        description = line[5].rstrip(".")
        entry["tool"] = [line[3]]
        entry["start"] = int(line[6])
        entry["end"] = int(line[7])
        
        # normalise zinc finger domain types, since these can have many alternative
        # descriptors
        if "zinc finger" in description.lower():
            description = "Zinc finger"
        
        if len(line) == 11:
            entry["domain_type"] = set([description])
        elif len(line) > 11:
            entry["domain_type"] = set([description, line[12]])
            entry["accession"] = line[11]
        
        # normalise the domain type from the "Coils" tool
        if entry["tool"] == ["Coils"]:
            entry["domain_type"] = set(["Coiled coil"])
        
        entry["domain_type"].discard("")
        
        return entry
    
    def parse_interpro_results(self, api_response):
        """ parses the output of an InterPro analysis run
        """
        
        domains = []
        for line in api_response.strip().split("\n"):
            line = line.split("\t")
            
            entry = self.build_domain_entry(line)
            
            # I'm not using PANTHER, Gene3D or SUPERFAMILY results, since Ensembl
            # doesn't include them (and because they have slightly different naming
            # conventions, whichmakes more difficult to match predictions between
            # different tools)
            if len(set(entry["tool"]) & set(["PANTHER", "Gene3D", "SUPERFAMILY",
                    "SignalP_GRAM_POSITIVE", "SignalP_GRAM_NEGATIVE"])) > 0:
                continue
            
            match = self.find_domain_match(entry, domains)
            
            if match > -1:
                # We extend the original to contain the overlaps and tool info.
                domain = domains[match]
                domains[match]["start"] = min(entry["start"], domain["start"])
                domains[match]["end"] = max(entry["end"], domain["end"])
                domains[match]["domain_type"] = entry["domain_type"] | domain["domain_type"]
                domains[match]["tool"] = list(set(entry["tool"]) | set(domain["tool"]))
            else:
                domains.append(entry)
        
        return domains
