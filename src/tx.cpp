#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <stdexcept>
#include <iostream>

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
    if (!is_exonic(cds_min)) {
        Region region = fix_cds_boundary(cds_min);
        int delta = std::abs(region.end - region.start);
        
        cds[0] = Region {cds[0].start + delta, cds[0].end};
        cds_min = region.start;
        cds.insert(cds.begin(), region);
    }
    
    if (!is_exonic(cds_max)) {
        Region region = fix_cds_boundary(cds_max);
        int delta = std::abs(region.end - region.start);
        
        int idx = cds.size() - 1;
        cds[idx] = Region {cds[idx].start, cds[idx].end - delta};
        cds_max = region.end;
        cds.push_back(region);
    }
    
    cds_length = 0;
    for (auto region : cds) {
        cds_length += (region.end - region.start) + 1;
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
    
    if ( is_exonic(position) ) {
        throw std::invalid_argument( "position shouldn't be in exons" );
    }
    
    Region region = get_closest_exon(position);
    
    int start_dist = std::abs(position - region.start);
    int end_dist = std::abs(position - region.end);
    
    int start;
    int end;
    if (start_dist < end_dist) {
        // if we are closer to the start, then we go back an exon
        int i = closest_exon_num(region.start) - 1;
        end = exons[i].end;
        start = end - start_dist;
    } else {
        // if we are closer to the end, then we go forward an exon
        int i = closest_exon_num(region.end) + 1;
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

bool Tx::is_exonic(int pos) {
    /**
       checks if a position lies within the exon ranges
       
       @position integer chromosome position e.g. 10000000
    */
    
    Region exon = get_closest_exon(pos);
    return (pos >= exon.start) && (pos <= exon.end);
}

bool compareRegion(const Region& a, int b) {
    // include a sort operator for Region structs. This allows quick searching
    // of vectors of Regions
    return a.start < b;
}

int Tx::closest_exon_num(int pos) {
    /* find the index of the closest exon to a chromosome position
    */
    return closest_exon_num(pos, exons);
}

int Tx::closest_exon_num(int pos, std::vector<Region> group) {
    /**
        find the index for the closest exon/CDS to a chromosome position
        
        @position: chromosomal nucleotide position
        
        @returns number of exon containing the position
    */
    int idx = std::lower_bound(group.begin(), group.end(), pos, compareRegion) - group.begin();
    int size = group.size();
    if (idx == 0) {
        return idx;
    } else if (idx == size and pos > group[idx - 1].end) {
        return idx - 1;
    }
    
    // get the exon, and the preceding exon (lower bound finds exon after pos)
    Region a = group[idx - 1];
    Region b = group[idx];
    
    // if the position is within an exon, return that exon number
    if ((pos >= a.start) && (pos <= a.end)) {
        return idx - 1;
    } else if ((pos >= b.start) && (pos <= b.end)) {
        return idx;
    }
    
    // pos lies after the start of a. a is upstream, b is downstream, find closest.
    int delta_a = std::min(std::abs(pos - a.start), std::abs(pos - a.end));
    int delta_b = std::min(std::abs(pos - b.start), std::abs(pos - b.end));
    
    if (delta_a < delta_b) {
        return idx - 1;
    }
    
    return idx;
}

Region Tx::get_closest_exon(int pos) {
    /**
       find the closest exon for a chromosome position
       
       @position integer chromosome position e.g. 10000000
    */
    int idx = closest_exon_num(pos);
    return exons[idx];
}

bool Tx::in_coding_region(int pos) {
    /**
       checks if a position lies within the CDS
       
       @position integer chromosome position e.g. 10000000
    */
    int idx = closest_exon_num(pos, cds);
    Region region = cds[idx];
    return (pos >= region.start) && (pos <= region.end);
}

CDS_coords Tx::to_closest_exon(int pos) {
    /* shift an intronic site to the nearest exon boundary (but track the offset)
    */
    Region exon = get_closest_exon(pos);
    int site;
    if (std::abs(exon.start - pos) < std::abs(exon.end - pos)) {
        site = exon.start;
    } else {
        site = exon.end;
    }
    
    bool fwd = get_strand() == '+';
    int offset = (fwd) ? pos - site : site - pos;
    return CDS_coords {site, offset};
}

CDS_coords Tx::get_coding_distance(int pos) {
    /**
        find the distance in coding bp to the CDS start
        
        @pos_1: first chromosome position
        
        @returns CDS_coords, which contains distance in coding sequence base
            pairs between two positions and an offset for distance into the
            intron
    */
    int cds_start = get_cds_start();
    
    int offset = 0;
    if ( !is_exonic(pos) ) {
        CDS_coords site = to_closest_exon(pos);
        offset = site.offset;
        pos = site.position;
    }
    
    int first = closest_exon_num(cds_start);
    int alt = closest_exon_num(pos);
    
    bool fwd = get_strand() == '+';
    bool after = alt > first;
    bool positive = (fwd & after) | (!fwd & !after);
    // use positive coding position if the transcript is either:
    //   a) on the + strand and the site is 3' of CDS start
    //   b) on the - strand and the site is 5' of CDS start
    int delta = 0;
    if (first == alt) {
        delta = (fwd) ? pos - cds_start : cds_start - pos ;
    } else {
        int dist_1;
        int dist_2;
        if (after) {
            dist_1 = (fwd) ? exons[first].end - cds_start : cds_start - exons[first].end;
            dist_2 = (fwd) ? (pos - exons[alt].start) + 1 : (exons[alt].start - pos) + 1;
        } else {
            dist_1 = (fwd) ? exons[first].start - cds_start : cds_start - exons[first].start;
            dist_2 = (fwd) ? (pos - exons[alt].end) - 1 : (exons[alt].end - pos) + 1;
        }
        
        int lo = std::min(first, alt) + 1;
        int hi = std::max(first, alt);
        delta = dist_1 + dist_2;
        for (int i=lo; i < hi; i++) {
            int left = exons[i].start;
            int right = exons[i].end;
            delta += (positive) ? (right - left) + 1 : (left - right) - 1;
        }
    }
    
    return CDS_coords {delta, offset};
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
        // get exon boundaries as CDS distance from the start site
        CDS_coords start = get_coding_distance(region.start);
        CDS_coords end = get_coding_distance(region.end);
        
        // cache the CDS positions of the exon boundaries
        exon_to_cds[region.start] = std::abs(start.position);
        exon_to_cds[region.end] = std::abs(end.position);
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
                return start + (end_cds - cds_position) - offset;
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
    
    cds_length = 0;
    for (auto region : cds) {
        cds_length += (region.end - region.start) + 1;
    }
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
    if (length % 2 != 1) {
        throw std::invalid_argument( "length is not an odd number") ;
    }
    
    int start = pos - floor(length/2);
    int end = start + length;
    return get_seq_in_region(start, end);
}

std::string Tx::get_seq_in_region(int start, int end) {
    if (start < get_start() - gdna_offset || end > get_end() + gdna_offset) {
        throw std::invalid_argument( "sequence position not in gene range" );
    }
    int pos = start - get_start() + gdna_offset;
    int length = end - start;
    std::string seq = genomic_sequence.substr(pos, length);
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    return seq;
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
    CDS_coords site = get_coding_distance(bp);
    
    if ((site.position < 0) | (site.position > cds_length)) {
        throw std::invalid_argument( "position not inside CDS region" );
    }
    
    // define the default values for the codon positions. We check these later
    // in python and convert to None if they are still the defaults.
    int codon_number = -9999999;
    int intra_codon = -1;
    std::string codon_seq = "";
    std::string initial_aa = "";
    
    CDS_coords cds_pos = get_coding_distance(bp);
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
    
    Region exon = get_closest_exon(bp);
    int distance = std::min(std::abs(exon.start - bp), std::abs(exon.end - bp));
    
    // sites within the coding region are actually one bp further away,
    // since we are measuring the distance to the base inside the exon
    // boundary
    if ( in_coding_region(bp) ) {
        distance += 1;
    }
    
    return distance;
}

bool Tx::overlaps_cds(int start, int end) {
    /* checks if a genomic range overlaps any exon range
    */
    
    // find the
    int idx = closest_exon_num(start, cds);
    int prev = std::max(0, idx - 1);
    int post = std::min((int) cds.size() - 1, idx + 1);
    
    if (start <= cds[idx].end && cds[idx].start < end) {
        return true;
    } else if (start < cds[prev].end && cds[prev].start < end) {
        return true;
    } else if (start < cds[post].end && cds[post].start < end) {
        return true;
    }
    return false;
}

std::string Tx::outside_gene_cq(int start, int end, std::string alt) {
    bool fwd = get_strand() == '+';
    if ((end < tx_start) & ((tx_start - end) < 5000)) {
        if (fwd) {
            return "upstream_gene_variant";
        } else {
            return "downstream_gene_variant";
        }
    } else if ((start > tx_end) & ((start - tx_end) < 5000)) {
        if (fwd) {
            return "downstream_gene_variant";
        } else {
            return "upstream_gene_variant";
        }
    }
    return "intergenic_variant";
}

std::string Tx::intronic_cq(int start, int end) {
    /* figure out consequence for a site not in coding region
    */
    // figure out if we are in a very short intron
    auto idx = closest_exon_num(start);
    auto exon = exons[idx];
    int intron_length;
    if (exon.end < start) {
        auto other = exons[idx + 1];
        intron_length = other.start - exon.end;
    } else {
        auto other = exons[idx - 1];
        intron_length = exon.start - other.end;
    }
    if (intron_length < 10) {
        return "coding_sequence_variant";
    }
    
    // TODO: figure out for variants spanning an exon
    CDS_coords start_coord = get_coding_distance(start);
    CDS_coords end_coord = get_coding_distance(end);
    
    bool intronic = (start_coord.offset != 0) & (end_coord.offset != 0);
    if (!intronic) {
        if (end_coord.position <= 0) {
            return "5_prime_UTR_variant";
        } else if (start_coord.position >= cds_length) {
            return "3_prime_UTR_variant";
        }
    }
    
    // ok, now we have intronic variants, so we check for splice variants first
    int delta = std::min(std::abs(start_coord.offset), std::abs(end_coord.offset));
    if ((delta < 3) & (start_coord.offset < 0)) {
        return "splice_acceptor_variant";
    } else if ((delta < 3) & (start_coord.offset > 0)) {
        return "splice_donor_variant";
    } else if (delta < 9) {
        return "splice_region_variant";
    }
    
    return "intron_variant";
}

std::string Tx::coding_cq(int start, int end, std::string alt) {
    if (get_strand() == '-') {
        alt = reverse_complement(alt);
    }
    // TODO: figure out initial and mutated amino acids
    Codon codon = get_codon_info(start);
    std::string initial_aa = codon.initial_aa;
    codon.codon_seq.replace(codon.intra_codon, 1, alt);
    std::string mutated_aa = translate(codon.codon_seq);
    
    if (initial_aa != "*" && mutated_aa == "*") { return "stop_gained"; }
    if (initial_aa == "*" && mutated_aa != "*") { return "stop_lost"; }
    
    if (initial_aa != mutated_aa) {
        if (codon.codon_number == 0) {
            return "start_lost";
        }
        return "missense_variant";
    }
    
    if (codon.codon_number == 0) {
        return "start_retained_variant";
    } else if (get_codon_number_for_cds_position(cds_length - 1) == codon.codon_number) {
        return "stop_retained_variant";
    }
    return "synonymous_variant";
}

int min_len(std::string a, std::string b) {
    return std::min(a.size(), b.size());
}

void Tx::trim_alleles(int& start, int& end, std::string& alt) {
    /* realign alt against ref sequence, and strip down to alternate bases
    */
    std::string ref = get_seq_in_region(start, end + 1);
    // trim the front of the alleles
    while ((min_len(ref, alt) > 0) && (ref.size() > 1) && (alt[0] == ref[0])) {
        ref.erase(0, 1);
        alt.erase(0, 1);
        start += 1;
    }
    
    // trim the back of the alleles
    while ((min_len(ref, alt) > 0) && (ref.size() > 1) && (alt.back() == ref.back())) {
        ref.pop_back();
        alt.pop_back();
        end -= 1;
    }
}

std::string Tx::consequence(int start, int end, std::string alt) {
    if (start > end) {
        throw std::invalid_argument("start position is less than end position");
    }
    
    if ((end <= tx_start) | (start >= tx_end)) {
        return outside_gene_cq(start, end, alt);
    }
    
    // TODO: need to realign alt against ref sequence, and strip down to alternate bases
    trim_alleles(start, end, alt);
    
    if (!overlaps_cds(start, end)) {
        return intronic_cq(start, end);
    }
    
    // at this point we know that the site lies within the gene boundary, and
    // at least partially overlaps the CDS
    std::string ref = get_seq_in_region(start, end+1);
    int ref_len = ref.size();
    int alt_len = alt.size();
    int delta = ref_len - alt_len;
    if (ref.size() != alt.size()) {
        if ((delta % 3) != 0) {
            return "frameshift_variant";
        } else if (ref_len > alt_len) {
            return "inframe_deletion";
        } else if (alt_len > ref_len) {
            return "inframe_insertion";
        }
    }
    
    return coding_cq(start, end, alt);
}
