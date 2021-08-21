
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "gzstream/gzstream.h"

#include "gtf.h"

namespace gencode {

// parse the required firleds from the attributes field
void get_attributes_fields(GTFLine &info, std::string &line, int offset) {
    // we could check for each field individually, but since we know the order
    // of the fields, it's much quicker to just search the remaining substring
    int tx_start = line.find("transcript_id", offset) + 15;
    int tx_end = line.find("\"", tx_start);

    if (tx_start - 15 == std::string::npos) {
        tx_start = tx_end;  // handle if the string was not found
    }

    int gene_start = line.find("gene_name", tx_end) + 11;
    int gene_end = line.find("\"", gene_start);

    int type_start = line.find("transcript_type", gene_end) + 17;
    int type_end = line.find("\"", type_start);

    if (type_start - 17 == std::string::npos) {
        type_start = type_end;  // handle if the string was not found
    }
    
    info.symbol = line.substr(gene_start, gene_end - gene_start);
    info.tx_id = line.substr(tx_start, tx_end - tx_start);
    info.transcript_type = line.substr(type_start, type_end - type_start);
}

// parse required fields from a GTF line
GTFLine parse_gtfline(std::string & line) {
    if (line.size() == 0) {
        throw std::out_of_range("end of file");
    }

    GTFLine info;

    // there are only a few fields we need from the GTF lines, and some fields
    // are only a single character long, so it's quickest to search for the next
    // tab along the line, then extract the substring to get the required fields.
    // getline() with tab delimiter was 2X slower.
    int chr_idx = 0;
    int source_idx = line.find("\t", chr_idx + 4);
    int feature_idx = line.find("\t", source_idx + 6);
    int start_idx = line.find("\t", feature_idx + 3);
    int end_idx = line.find("\t", start_idx + 3);
    int score_idx = line.find("\t", end_idx + (end_idx - start_idx));

    info.chrom = line.substr(chr_idx, source_idx - chr_idx);
    info.feature = line.substr(feature_idx + 1, start_idx - feature_idx - 1);
    info.start = std::stoi(line.substr(start_idx + 1, end_idx - start_idx - 1));
    info.end = std::stoi(line.substr(end_idx + 1, score_idx - end_idx - 1));
    info.strand = line[score_idx + 3];

    get_attributes_fields(info, line, score_idx + 6);

    return info;
}

// open GTF file handle
GTF::GTF(std::string path) {
    gzipped = path.substr(path.length()-2, 2) == "gz";
    if (gzipped) {
        gzhandle.open(path.c_str());
    } else {
        handle.open(path, std::ios::in);
    }
}

// get next line from the GTF
GTFLine GTF::next() {
    if (gzipped) { 
        std::getline(gzhandle, line);
    } else { 
        std::getline(handle, line);
    }
    while (line[0] == '#') {
        if (gzipped) {
            std::getline(gzhandle, line);
        } else {
            std::getline(handle, line);
        }
    }
    return parse_gtfline(line);
}

} // namespace
