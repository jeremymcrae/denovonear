#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "tx.h"
#include "weighted_choice.h"

Region get_gene_range(Tx tx) {
    /**
        get the lowest and highest positions of a transcripts coding sequence
    */
    
    int boundary_1 = tx.get_cds_start();
    int boundary_2 = tx.get_cds_end();
    
    int start = min(boundary_1, boundary_2);
    int end = max(boundary_1, boundary_2);
    
    return Region {start, end}
}

std::string get_mutated_aa(Tx tx, std::string base, std::string codon, int intra_codon) {
    /**
        find the amino acid resulting from a base change to a codon
        
        @tx transcript object for a gene
        @base alternate base (e.g. 'G') to introduce
        @codon DNA sequence of a single codon
        @intra_codon position within the codon to be altered (0-based)
    
        @returns single character amino acid code translated from the altered
            codon.
    */
      
    // figure out what the mutated codon is
    codon.replace(intra_codon, 1, base);
    
    return tx.translate(codon);
}

void SitesChecks::SitesChecks(Tx tx, std::vector<std::vector<std::string>> rates,
    Tx masked_sites) {
    
    _tx = tx;
    masked = masked_sites;
    
    // convert the rates data to a nested map
    int len = rates.size();
    for (int i=1; i < len, i++) {
        std::vector<std::string> line = rates[i];
        rates[line[0]][line[1]] = std::stod(line[2]);
    }
    
    // initialise a WeightedChoice object for each consequence category
    len = categories.size();
    std::map<std::string, WeightedChoice> rates;
    for (int i=1; i < len, i++) {
        std::string cq = categories[i];
        rates[cq] = WeightedChoice();
    }
    
    // check the consequence alternates for each base in the coding sequence
    Region region = get_gene_range(tx);
    for (int i=region.start; i < region.end + 1; i++ ) {
        check_position(i)
    }
}

bool SitesChecks::splice_lof_check(std::string initial_aa, std::string mutated_aa, int position) {
    /**
        checks if a variant has a splice_donor or splice_acceptor consequence
        
        These variants are defined as being the two intronic base-pairs adjacent
        to the intron/exon boundary.
    */
    
    return (!_tx.in_coding_region(position)) && boundary_dist < 3;
}

bool SitesChecks::nonsense_check(std::string initial_aa, std::string mutated_aa, int position) {
    /**
        checks if two amino acids are a nonsense (eg stop_gained) mutation
    */
    
    return initial_aa != "*" & mutated_aa == "*";
}

bool SitesChecks::missense_check(std::string initial_aa, std::string mutated_aa, int position) {
    /**
        checks if two amino acids are a missense mutation (but not nonsense)
    */
    
    // trim out nonsense mutations such as stop_gained mutations, and splice
    // site mutations
    if (nonsense_check(initial_aa, mutated_aa, position) || \
            splice_lof_check(initial_aa, mutated_aa, position)) {
        return false;
    }
    
    // include the site if it mutates to a different amino acid.
    return initial_aa != mutated_aa;
}

bool SitesChecks::splice_region_check(std::string initial_aa, std::string mutated_aa, int position) {
    /**
        checks if a variant has a splice_region consequence, but not
        splice_donor or splice_acceptor
    */
    
    if ( splice_lof_check(initial_aa, mutated_aa, position) || \
            initial_aa != mutated_aa ) {
        return false;
    }
    
    // catch splice region variants within the exon, and in the appropriate
    // region of the intron (note that loss of function splice_donor and
    // splice_acceptor variants have been excluded when we trimmed nonsense).
    if ( _tx.in_coding_region(position) ) {
        // check for splice_region_variant inside exon
        return boundary_dist < 4;
    } else {
        // check for splice_region_variant inside intron
        return boundary_dist < 9;
    }
}

bool SitesChecks::synonymous_check(std::string initial_aa, std::string mutated_aa, int position) {
    /**
        checks if two amino acids are synonymous
    */
    
    return !nonsense_check(initial_aa, mutated_aa, position)
        && !splice_lof_check(initial_aa, mutated_aa, position)
        && !missense_check(initial_aa, mutated_aa, position)
        && !splice_region_check(initial_aa, mutated_aa, position));
}

void SitesChecks::check_position(bp) {
    /**
        add the consequence specific rates for the alternates for a variant
        
        @bp genomic position of the variant
    */
    
    // ignore sites within masked regions (typically masked because the
    // site has been picked up on alternative transcript)
    if ( masked != None and masked.in_coding_region(bp) ) {
        return ;
    }
    
    // ignore sites outside the CDS region
    if (bp < min(_tx.get_cds_start(), _tx.get_cds_end()) ||
        bp > max(_tx.get_cds_start(), _tx.get_cds_end())) {
        return ;
    }
    
    std::string seq = _tx.get_trinucleotide(bp);
    boundary_dist = _tx.get_boundary_distance(bp);
    
    char fwd = "+";
    if (_tx.get_strand() == fwd) {
        seq = _tx.reverse_complement(seq);
    }
    
    try {
        Codon codon = _tx.get_codon_info(bp);
    } catch {
        return ;
    }
    
    std::string initial_aa = codon.initial_aa;
    int cds_pos = codon.cds_pos;
    
    // drop the initial base, since we want to mutate to other bases
    std::vector<std::string> alts = bases;
    bases.erase(std::find(seq[1]));
    
    for (int i=1; i < 4; i++) {
        std::string base = alts[i];
        std::string mutated_aa = initial_aa;
        std::string alt_base = seq[0] + base + seq[2];
        rate = mut_dict[seq][alt_base]
        if ( initial_aa != "" ) {
            std::string mutated_aa = get_mutated_aa(_tx, base, codon.codon_seq, codon.intra_codon)
        }
        
        std::string category = "";
        if (nonsense_check(initial_aa, mutated_aa, bp)) {
            category = "nonsense";
        } else if (splice_lof_check(initial_aa, mutated_aa, bp)) {
            category = "splice_lof";
        } else if (missense_check(initial_aa, mutated_aa, bp)) {
            category = "missense";
        } else if (splice_region_check(initial_aa, mutated_aa, bp)) {
            category = "splice_region";
        } else if (synonymous_check(initial_aa, mutated_aa, bp)) {
            category = "synonymous";
        }
        
        // figure out what the ref and alt alleles are, with respect to
        // the + strand.
        std::string ref_base = seq[1];
        std::string alt = base;
        if (_tx.get_strand() != fwd) {
            ref_base = transdict[ref_base];
            alt = transdict[alt];
        }
        
        rates[category].add_choice(cds_pos, rate, ref_base, alt);
        
        if (category == "nonsense" || category == "splice_lof"]){
            rates["loss_of_function"].add_choice(cds_pos, rate, ref_base, alt);
        }
    }
}
