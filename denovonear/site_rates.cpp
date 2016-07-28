#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

// #include "tx.h"
// #include "weighted_choice.h"
#include "site_rates.h"

Region _get_gene_range(Tx tx) {
    /**
        get the lowest and highest positions of a transcripts coding sequence
    */
    
    int boundary_1 = tx.get_cds_start();
    int boundary_2 = tx.get_cds_end();
    
    int start = std::min(boundary_1, boundary_2);
    int end = std::max(boundary_1, boundary_2);
    
    return Region {start, end};
}

std::string _get_mutated_aa(Tx tx, std::string base, std::string codon, int intra_codon) {
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

void SitesChecks::init(std::vector<std::vector<std::string>> mut) {
    
    // convert the rates data to a nested map
    int len = mut.size();
    for (int i=0; i < len; i++) {
        std::vector<std::string> line = mut[i];
        mut_dict[line[0]][line[1]] = std::stod(line[2]);
    }
    
    initialise_choices();
    // check the consequence alternates for each base in the coding sequence
    Region region = _get_gene_range(_tx);
    for (int i=region.start; i < region.end + 1; i++ ) {
        check_position(i);
    }
}

void SitesChecks::initialise_choices() {
    // initialise a WeightedChoice object for each consequence category
    for (auto it=categories.begin(); it != categories.end(); it++) {
        rates[*it] = Chooser();
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
    
    return (initial_aa != "*") & (mutated_aa == "*");
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
    return (initial_aa != mutated_aa);
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
    
    return !nonsense_check(initial_aa, mutated_aa, position) && \
         !splice_lof_check(initial_aa, mutated_aa, position) && \
         !missense_check(initial_aa, mutated_aa, position) && \
         !splice_region_check(initial_aa, mutated_aa, position);
}

void SitesChecks::check_position(int bp) {
    /**
        add the consequence specific rates for the alternates for a variant
        
        @bp genomic position of the variant
    */
    
    // // ignore sites within masked regions (typically masked because the
    // // site has been picked up on alternative transcript)
    // if ( has_mask && masked.in_coding_region(bp) ) {
    //     return ;
    // }
    
    // ignore sites outside the CDS region
    if (bp < std::min(_tx.get_cds_start(), _tx.get_cds_end()) ||
        bp > std::max(_tx.get_cds_start(), _tx.get_cds_end())) {
        return ;
    }
    
    std::string seq = _tx.get_trinucleotide(bp);
    boundary_dist = _tx.get_boundary_distance(bp);
    
    char fwd = '+';
    if (_tx.get_strand() != fwd) {
        seq = _tx.reverse_complement(seq);
    }
    
    Codon codon;
    try {
        codon = _tx.get_codon_info(bp);
    } catch ( const std::invalid_argument& e ) {
        return ;
    }
    
    std::string initial_aa = codon.initial_aa;
    int cds_pos = codon.cds_pos;
    
    // drop the initial base, since we want to mutate to other bases
    std::string ref = seq.substr(1, 1);
    std::vector<std::string> alts(bases);
    alts.erase(std::find(alts.begin(), alts.end(), ref));
    
    for (auto it=alts.begin(); it != alts.end(); it++) {
        std::string alt = *it;
        std::string mutated_aa = initial_aa;
        std::string alt_seq = seq[0] + alt + seq[2];
        double rate = mut_dict[seq][alt_seq];
        if ( initial_aa != "" ) {
            mutated_aa = _get_mutated_aa(_tx, alt, codon.codon_seq, codon.intra_codon);
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
        ref = seq.substr(1, 1);
        if (_tx.get_strand() != fwd) {
            ref = transdict[ref];
            alt = transdict[alt];
        }
        
        rates[category].add_choice(cds_pos, rate, ref, alt);
        
        if (category == "nonsense" || category == "splice_lof") {
            rates["loss_of_function"].add_choice(cds_pos, rate, ref, alt);
        }
    }
}
