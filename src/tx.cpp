#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <stdexcept>

#include "tx.h"

Tx::Tx(std::string transcript_id, std::string chromosome,
    int start_pos, int end_pos, char strand) {
    /**
        Constructor for Tx class
    */
    
    name = transcript_id;
    
    chrom = chromosome;
    tx_start = start_pos;
    tx_end = end_pos;
    
    if ( strand != '+' && strand != '-' ) {
        throw std::invalid_argument( "unknown strand type" );
    }
    tx_strand = strand;
}

void Tx::set_exons(std::vector<std::vector<int>> exon_ranges,
    std::vector<std::vector<int>> cds_ranges) {
    /**
       set exon ranges to the class object
       
       @exon_ranges list of lists e.g. [[5, 10], [20, 30]]
    */
    
    exons.clear();
    std::sort(exon_ranges.begin(), exon_ranges.end());
    
    for (auto range : exon_ranges) {
        Region region {range[0], range[1]};
        exons.push_back(region);
    }
    
    // If the transcript lacks exon coordinates and only has a single CDS
    // region, then if the CDS region fits within the gene range, make a
    // single exon, using the transcript start and end. This prevents issues
    // when adding the CDS coordinates for the transcript.
    if (exons.empty()) {
        std::vector<int> exon = cds_ranges[0];
        if (cds_ranges.size() == 1 && std::min(exon[0], exon[1]) >= tx_start \
                && std::max(exon[0], exon[1]) <= tx_end) {
            Region region {tx_start, tx_end};
            exons.push_back(region);
        } else {
            std::string msg = get_name() + " lacks exon coordinates";
            throw std::invalid_argument( msg );
        }
    }
}

void Tx::set_cds(std::vector<std::vector<int>> cds_ranges) {
    /**
       set CDS ranges to the class object
       
       @cds_ranges nested vector of ints e.g. [[5, 10], [20, 30]]
   */
   
   cds.clear();
   
    if ( exons.empty() ) {
        throw std::invalid_argument( "you need to set the exons before the CDS" );
    }
    
    std::sort(cds_ranges.begin(), cds_ranges.end());
    
    for (auto range : cds_ranges) {
        Region region {range[0], range[1]};
        cds.push_back(region);
    }
    
    cds_min = cds[0].start;
    cds_max = cds.back().end;
    
    // sometimes the cds start and end lie outside the exons. We find the
    // offset to the closest exon boundary, and use that to place the
    // position within the up or downstream exon, and correct the initial
    // cds boundaries
    if (!in_exons(cds_min)) {
        Region region = fix_cds_boundary(cds_min);
        int delta = std::abs(region.end - region.start);
        
        cds[0] = Region {cds[0].start + delta, cds[0].end};
        cds_min = region.start;
        cds.insert(cds.begin(), region);
    }
    
    if (!in_exons(cds_max)) {
        Region region = fix_cds_boundary(cds_max);
        int delta = std::abs(region.end - region.start);
        
        int idx = cds.size() - 1;
        cds[idx] = Region {cds[idx].start, cds[idx].end - delta};
        cds_max = region.end;
        cds.push_back(region);
    }
}

Region Tx::fix_cds_boundary(int position) {
    /**
        fix CDS start and ends which occur outside exons.
        
        The main cause is when the stop codon overlaps an exon boundary. In
        some of these cases, the CDS end is placed within the adjacent intron
        sequence (if that allows stop site sequence). To fix this, we adjust the
        position to the offset distance within the adjacent exon.
        
        @position chromosome position
        
        @returns adjusted chromosome positions
    */
    
    if ( in_exons(position) ) {
        throw std::invalid_argument( "position shouldn't be in exons" );
    }
    
    Region region = find_closest_exon(position);
    
    int start_dist = std::abs(position - region.start);
    int end_dist = std::abs(position - region.end);
    
    int start;
    int end;
    if (start_dist < end_dist) {
        // if we are closer to the start, then we go back an exon
        int i = get_exon_containing_position(region.start, exons) - 1;
        end = exons[i].end;
        start = end - start_dist;
    } else {
        // if we are closer to the end, then we go forward an exon
        int i = get_exon_containing_position(region.end, exons) + 1;
        start = exons[i].start;
        end = start + end_dist;
    }
    
    return Region {start, end};
}

int Tx::get_cds_start() {
    
    char fwd = '+';
    if (get_strand() == fwd) {
        return cds_min;
    } else {
        return cds_max;
    }
}

int Tx::get_cds_end() {
    
    char fwd = '+';
    if (get_strand() == fwd) {
        return cds_max;
    } else {
        return cds_min;
    }
}

bool Tx::in_exons(int position) {
    /**
       checks if a position lies within the exon ranges
       
       @position integer chromosome position e.g. 10000000
    */
    
    for (auto &region : exons) {
        if (region.start <= position && position <= region.end) {
            return true;
        }
    }
    
    return false;
}

Region Tx::find_closest_exon(int position) {
    /**
       checks if a position lies within the exon ranges
       
       @position integer chromosome position e.g. 10000000
    */
    
    return find_closest_exon(position, exons);
}

Region Tx::find_closest_exon(int position, std::vector<Region> & ranges) {
    /**
       checks if a position lies within the exon ranges
       
       @position integer chromosome position e.g. 10000000
    */
    
    long max_distance = 1000000000;
    int ref_start = 0;
    int ref_end = 0;
    
    for (auto &region : ranges) {
        int start_dist = std::abs(region.start - position);
        int end_dist = std::abs(region.end - position);
        
        if (start_dist < max_distance) {
            ref_start = region.start;
            ref_end = region.end;
            max_distance = start_dist;
        }
        
        if (end_dist < max_distance) {
            ref_start = region.start;
            ref_end = region.end;
            max_distance = end_dist;
        }
    }
    
    return Region {ref_start, ref_end};
}

bool Tx::in_coding_region(int position) {
    /**
       checks if a position lies within the exon ranges
       
       @position integer chromosome position e.g. 10000000
    */
    
    for (auto &region : cds) {
        if (region.start <= position && position <= region.end) {
            return true;
        }
    }
    
    return false;
}

int Tx::get_exon_containing_position(int position, std::vector<Region> & ranges) {
    /**
        find the exon number for a position within an exon
        
        @position: chromosomal nucleotide position
        @ranges: vector of exon regions
        
        @returns number of exon containing the position
    */
    
    int exon_num = 0;
    
    for (auto &region : ranges) {
        if (region.start <= position && position <= region.end) {
            return exon_num;
        }
        exon_num += 1;
    }
    
    throw std::logic_error( "you've tried to identify the region containing a"
        "position that doesn't occur within the defined set of regions" );
}

int Tx::get_coding_distance(int pos_1, int pos_2) {
    /**
        find the distance in coding bp between two chr positions
        
        @pos_1: first chromosome position
        @pos_2: second chromosome position
        
        @returns distance in coding sequence base pairs between two positions
        
        @raises invalid_argument if position not in coding region
    */
    
    // make sure that the positions are within the coding region, otherwise
    // there's no point trying to calculate the coding distance
    if ( !in_coding_region(pos_1) ) {
        throw std::invalid_argument( "not in coding region" );
    }
    if ( !in_coding_region(pos_2) ) {
        throw std::invalid_argument( "not in coding region" );
    }
    
    int min_pos = std::min(pos_1, pos_2);
    int max_pos = std::max(pos_1, pos_2);
    
    int exon_1 = get_exon_containing_position(min_pos, cds);
    int exon_2 = get_exon_containing_position(max_pos, cds);
    
    int delta = 0;
    if (exon_1 == exon_2) {
        delta = max_pos - min_pos;
    } else {
        int dist_1 = cds[exon_1].end - min_pos;
        int dist_2 = (max_pos - cds[exon_2].start) + 1;
        
        delta = dist_1 + dist_2;
        for (int i=exon_1 + 1; i < exon_2; i++) {
            delta += (cds[i].end - cds[i].start) + 1;
        }
    }
    
    return delta;
}

CDS_coords Tx::chrom_pos_to_cds(int pos) {
    /**
        returns a chromosome position as distance from CDS ATG start
    */
    
    // need to convert the de novo event positions into CDS positions
    int cds_start = get_cds_start();
    
    try {
        return CDS_coords { get_coding_distance(cds_start, pos), 0 } ;
    } catch ( const std::invalid_argument& e ) {
        // catch the splice site functional mutations
        Region exon = find_closest_exon(pos);
        
        int site;
        if (std::abs(exon.start - pos) < std::abs(exon.end - pos)) {
            site = exon.start;
        } else {
            site = exon.end;
        }
        
        int offset = pos - site;
        
        // catch variants near an exon, but where the exon isn't in the CDS
        if (!in_coding_region(site)) {
            std::string msg = "Not near coding exon: " + std::to_string(pos) +
                " in transcript " + get_name();
            throw std::logic_error(msg);
        }
        
        // ignore positions outside the exons that are too distant from a boundary
        if (std::abs(offset) >= 9) {
            std::string msg = "distance to exon (" + std::to_string(std::abs(offset)) +
                ") > 8 bp for " + std::to_string(pos) + " in transcript "+ get_name();
            throw std::logic_error(msg);
        }
        
        return CDS_coords { get_coding_distance(cds_start, site), offset };
    }
}

void Tx::_cache_exon_cds_positions() {
    /**
        cache the exon boundary positions as CDS coordinates.
        
        Rather than recalculate the CDS positions each time we want to convert a
        chromosome position to a CDS position, we calculate them once and cache
        the coordinates in a dictionary for rapid access.
    */
    
    exon_to_cds.clear();
    
    for (auto &region : cds) {
        // get the positions of the exon boundaries in CDS distance from
        // the start site
        int start_cds = get_coding_distance(get_cds_start(), region.start);
        int end_cds = get_coding_distance(get_cds_start(), region.end);
        
        // cache the CDS positions of the exon boundaries
        exon_to_cds[region.start] = start_cds;
        exon_to_cds[region.end] = end_cds;
    }
}

int Tx::get_position_on_chrom(int cds_position, int offset) {
    /**
        figure out the chromosome position of a CDS site
    
        @cds_position position of a variant in CDS.
        @returns chromosome bp position of the CDS site
    */
    
    // cache the exon boundaries in CDS distances
    if (exon_to_cds.empty()) { _cache_exon_cds_positions(); }
    
    int start;
    int end;
    int start_cds;
    int end_cds;
    
    // quickly find the exon containing the CDS position
    for (auto &region : cds) {
        start = region.start;
        end = region.end;
        
        start_cds = std::min(exon_to_cds[start], exon_to_cds[end]);
        end_cds = std::max(exon_to_cds[start], exon_to_cds[end]);
        
        if (start_cds <= cds_position && cds_position <= end_cds) {
            // convert the CDS position to a chromosomal coordinate
            char fwd = '+';
            if (get_strand() == fwd) {
                return start + (cds_position - start_cds) + offset;
            } else {
                return start + (end_cds - cds_position) + offset;
            }
        }
    }
    
    throw std::invalid_argument( "position not in CDS" );
}

int Tx::get_codon_number_for_cds_position(int cds_position) {
    /**
        figure out the codon position for a position
    */
    
    return floor((float)cds_position / 3);
}

int Tx::get_position_within_codon(int cds_position) {
    /**
        get the position within a codon (in 0 based format eg 0, 1, 2)
    */
    
    return cds_position % 3;
}

void Tx::add_cds_sequence(std::string cds_dna) {
    /**
        add the CDS sequence
    */
    
    cds_sequence = cds_dna;
}

void Tx::add_genomic_sequence(std::string gdna, int offset=0) {
    /**
        add and process the genomic sequence into a CDS sequence
        
        @gdna string for genomic sequence. If the transcript strand is '-',
            then the DNA sequence will be for the - strand, so we need to
            reorient the DNA to the + strand.
    */
    
    if ( cds.size() == 0 ) {
        std::string msg = "You need to add CDS coordinates before adding"
            "genomic sequence!";
        throw std::invalid_argument(msg);
    }
    
    char fwd = '+';
    if (get_strand() != fwd) {
        // orient the DNA sequence to the + strand.
        gdna = reverse_complement(gdna);
    }
    
    gdna_offset = offset;
    genomic_sequence = gdna;
    
    std::string cds_seq;
    for (auto &region : cds) {
        int x = std::abs(region.start - get_start()) + offset;
        int len = (region.end - region.start) + 1;
        std::string bases = genomic_sequence.substr(x, len);
        
        if (get_strand() == fwd) {
            cds_seq += bases;
        } else {
            cds_seq = reverse_complement(bases) + cds_seq;
        }
    }
    
    // don't check the CDS matches expectations if the CDS sequence isn't set
    if (cds_sequence == "") { cds_sequence = cds_seq; return ; }
    
    // do a sanity check to check that we've got the right cds sequence, this
    // fails for at least one gene (CCDC18), which begins with a N, and
    // throws the coordinates off
    if (cds_seq != cds_sequence) {
        std::string msg = "Coding sequence from gene coordinates doesn't match "
            "coding sequence obtained from Ensembl.\nTranscript:" + get_name() +
            "\n" + cds_seq + "\n\nshould be\n" + cds_sequence + "\n";
        throw std::invalid_argument(msg);
    }
    
    _fix_cds_length();
}

void Tx::_fix_cds_length() {
    /**
        correct the coding sequence of a transcript, if it misses bases.
        
        Some transcripts don't cover the full length of a gene, they terminate
        in the middle of an amino acid. This is rare, and it is rarer that we
        select these transcripts for analysis. TNS3 is an example, where the
        de novos in the gene can only be contained within a transcript that is
        incomplete at the 3' end. This causes problems when we try to translate
        the final codon.
        
        We simply extend the coding sequence 1-2 bases to make a complete codon.
    */
    
    int diff = cds_sequence.size() % 3;
    int end = get_cds_end();
    
    int last = cds.size() - 1;
    
    if (diff != 0) {
        diff = 3 - diff;
        
        char fwd = '+';
        if (get_strand() == fwd) {
            cds[last] = Region {cds[last].start, cds[last].end + diff};
            int start_bp = std::abs(end - get_start()) + gdna_offset;
            cds_sequence += genomic_sequence.substr(start_bp, diff);
        } else {
            cds[0] = Region {cds[0].start - diff, cds[0].end};
            int start_bp = std::abs(get_end() - end) + gdna_offset;
            cds_sequence += genomic_sequence.substr(start_bp, diff);
        }
    }
    
    cds_min = cds[0].start;
    cds_max = cds[last].end;
    
    // shifting the CDS coordinates can, very infrequently, shift the CDS
    // beyond the exon coordinates. This is only a problem if we take the union
    // transcript coordinates, but should be accounted for.
    if (cds_min < exons[0].start) {
        exons[0] = Region {exons[0].start - diff, exons[0].end};
    }
    
    last = exons.size() - 1;
    if (cds_max > exons[last].end) {
        exons[last] = Region {exons[last].start, exons[last].end + diff};
    }
}

std::string Tx::reverse_complement(std::string seq) {
    /**
        reverse complement a DNA or RNA sequence
    */
    
    std::reverse(seq.begin(), seq.end());
    std::string complement;
    
    for (auto &base : seq) {
        complement += transdict[base];
    }
    
    return complement;
}

std::string Tx::get_centered_sequence(int pos, int length) {
    /**
        obtains the sequence around a position
        
        @pos integer for chromosomal nucleotide position.
        @length integer for length of region to select.
        
        @returns nucleotide sequence e.g. a trinucleotide such as 'TCC',
            centered around the required position, given with respect to the fwd
            strand.
    */
    
    if (pos < 0) {
        throw std::invalid_argument( "sequence position < 0" );
    }
    if (pos <= get_start() - gdna_offset || pos >= get_end() + gdna_offset) {
        throw std::invalid_argument( "sequence position not in gene range" );
    }
    if (length % 2 != 1) {
        throw std::invalid_argument( "length is not an odd number") ;
    }
    
    int sequence_pos = pos - get_start() + gdna_offset - floor(length/2);
    
    return genomic_sequence.substr(sequence_pos, length);
}

std::string Tx::get_codon_sequence(int codon) {
    /**
        get the codon sequence for a given codon_number
    */
    
    if (codon < 0) {
        throw std::invalid_argument( "codon position < 0" );
    }
    if (codon > (int)cds_sequence.size() / 3) {
        throw std::invalid_argument( "codon position not in gene range" );
    }
    
    return cds_sequence.substr(codon * 3, 3);
}

std::string Tx::translate(std::string seq) {
    /**
        translate a DNA codon to a single character amino acid
        
        @seq codon sequence, or longer DNA sequences.
        
        @returns translated amino acid sequence as single character codes.
    */
    
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    
    std::string protein;
    int n = 3;
    int len = seq.size();
    std::string codon;
    
    for ( int i=0; i < len; i=i+n ) {
        codon = seq.substr(i, n);
        
        if (aa_code.count(codon) == 0) {
            std::string msg = "cannot translate codon: " + codon;
            throw std::invalid_argument( msg );
        }
        
        protein += aa_code[codon];
    }
    
    return protein;
}

Codon Tx::get_codon_info(int bp) {
    /**
        get the details of the codon which a variant resides in
        
        @bp nucleotide position of the variant (within the transcript gene range)
    
        @returns dictionary of codon sequence, cds position, amino acid that the
            codon translates to, and position within the codon.
    */
    
    bool in_coding = in_coding_region(bp);
    
    // ignore positions outside the exons that are too distant from a boundary
    if (not in_coding and get_boundary_distance(bp) >= 9) {
        throw std::invalid_argument( "position too distant from an exon" );
    }
    
    // define the default values for the codon positions. We check these later
    // in python and convert to None if they are still the defaults.
    int codon_number = -1;
    int intra_codon = -1;
    std::string codon_seq = "";
    std::string initial_aa = "";
    
    CDS_coords cds_pos = chrom_pos_to_cds(bp);
    if (in_coding) {
        codon_number = get_codon_number_for_cds_position(cds_pos.position);
        intra_codon = get_position_within_codon(cds_pos.position);
        codon_seq = get_codon_sequence(codon_number);
        initial_aa = translate(codon_seq);
    }
    
    return Codon {cds_pos.position, codon_seq, intra_codon, codon_number,
        initial_aa, cds_pos.offset};
}

int Tx::get_boundary_distance(int bp) {
    /**
        get the distance in bp for a variant to the nearest exon boundary
        
        @bp nucleotide position of the variant (within the transcript gene range)
        @returns distance in base-pairs to the nearest exon boundary.
    */
    
    Region exon = find_closest_exon(bp);
    int distance = std::min(std::abs(exon.start - bp), std::abs(exon.end - bp));
    
    // sites within the coding region are actually one bp further away,
    // since we are measuring the distance to the base inside the exon
    // boundary
    if ( in_coding_region(bp) ) {
        distance += 1;
    }
    
    return distance;
}
