#include <string>
#include <vector>
#include <cmath>
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
    
    int len = exon_ranges.size();
    std::sort(exon_ranges.begin(), exon_ranges.end());
    
    for (int i=0; i < len; i++) {
        std::vector<int> range = exon_ranges[i];
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
            std::string msg = "%s lacks exon coordinates", get_name();
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
    
    int len = cds_ranges.size();
    std::sort(cds_ranges.begin(), cds_ranges.end());
    
    for (int i=0; i < len; i++) {
        std::vector<int> range = cds_ranges[i];
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
        int delta = abs(region.end - region.start);
        
        cds[0] = Region {cds[0].start + delta, cds[0].end};
        cds_min = region.start;
        cds.insert(cds.begin(), region);
    }
    
    if (!in_exons(cds_max)) {
        Region region = fix_cds_boundary(cds_max);
        int delta = abs(region.end - region.start);
        
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
    
    int start_dist = abs(position - region.start);
    int end_dist = abs(position - region.end);
    
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
    
    char x = '+';
    if (get_strand() == x) {
        return cds_min;
    } else {
        return cds_max;
    }
}

int Tx::get_cds_end() {
    
    char x = '+';
    if (get_strand() == x) {
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
    
    int len = exons.size();
    
    for (int i=0; i < len; i++) {
        Region region = exons[i];
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

Region Tx::find_closest_exon(int position, std::vector<Region> ranges) {
    /**
       checks if a position lies within the exon ranges
       
       @position integer chromosome position e.g. 10000000
    */
    
    int len = ranges.size();
    
    long max_distance = 1000000000;
    int ref_start = 0;
    int ref_end = 0;
    
    for (int i=0; i < len; i++) {
        Region region = ranges[i];
        int start_dist = abs(region.start - position);
        int end_dist = abs(region.end - position);
        
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
    
    int len = cds.size();
    
    for (int i=0; i < len; i++) {
        Region region = cds[i];
        if (region.start <= position && position <= region.end) {
            return true;
        }
    }
    
    return false;
}

int Tx::get_exon_containing_position(int position, std::vector<Region> ranges) {
    /**
        find the exon number for a position within an exon
        
        @position: chromosomal nucleotide position
        @ranges: vector of exon regions
        
        @returns number of exon containing the position
    */
    
    int len = ranges.size();
    int exon_num = 0;
    
    for (int i=0; i < len; i++) {
        Region region = ranges[i];
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

int Tx::chrom_pos_to_cds(int pos) {
    /**
        returns a chromosome position as distance from CDS ATG start
    */
    
    // need to convert the de novo event positions into CDS positions
    int cds_start = get_cds_start();
    
    try {
        return get_coding_distance(cds_start, pos);
    } catch ( const std::invalid_argument& e ) {
        // catch the splice site functional mutations
        Region exon = find_closest_exon(pos);
        
        int site;
        if (abs(exon.start - pos) < abs(exon.end - pos)) {
            site = exon.start;
        } else {
            site = exon.end;
        }
        
        int distance = abs(site - pos);
        
        // catch variants near an exon, but where the exon isn't in the CDS
        if (!in_coding_region(site)) {
            std::string msg = "Not near coding exon: %i in transcript %s", pos, get_name();
            throw std::logic_error(msg);
        }
        
        // ignore positions outside the exons that are too distant from a boundary
        if (distance >= 9) {
            std::string msg = "distance to exon (%i) > 8 bp for %i in transcript %s", distance, pos, get_name();
            throw std::logic_error(msg);
        }
        
        return get_coding_distance(cds_start, site);
    }
}
