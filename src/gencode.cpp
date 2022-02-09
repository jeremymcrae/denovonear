
#include <algorithm>
#include <cstdint>
#include <string>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

#include <iostream>

#include "gtf.h"
#include "tx.h"
#include "gencode.h"

namespace gencode {

// check which exon is first, by start position
bool compareExons(std::vector<int> e1, std::vector<int> e2) {
    return (e1[0] < e2[0]);
}

void sort_exons(std::vector<std::vector<int> > & exons) {
    std::sort(exons.begin(), exons.end(), compareExons);
}

// find the index of the exon containing a given chromosome position
uint get_exon_num(std::vector<std::vector<int> > exons, int pos) {
    for (uint i=0; i<exons.size(); i++) {
        if (pos >= exons[i][0] && pos <= exons[i][1]) {
            return i;
        }
    }
    throw std::invalid_argument("can't find exon for an end codon");
}

// adjust CDS coords for any end codon positions
// 
// the CDS might exclude the start codon and end codon from the CDS range
// we have collected the start and end codon ends in cds_range. We can't
// just set the first CDS coord and last CDS coord to their values though,
// as at least one stop codon spans an intron boundary, which messes up the
// CDS if included as is.
void include_end_codons(std::map<std::string, int> cds_range, TxInfo & info) {
    if (info.cds.size() == 0) {
        return;
    }
    sort_exons(info.cds);
    sort_exons(info.exons);

    // handle left (5') boundary
    uint first_idx = get_exon_num(info.exons, info.cds[0][0]);
    uint min_idx = get_exon_num(info.exons, cds_range["min"]);
    if (min_idx == first_idx) {
        info.cds[0][0] = cds_range["min"];
    } else {
        info.cds[0][0] = info.exons[first_idx][0];  // extend existing CDS
        std::vector<int> extra_cds = {cds_range["min"], info.exons[min_idx][1]};
        auto it = info.cds.begin();
        info.cds.insert(it, extra_cds);
    }

    // handle right (3') boundary
    uint last_idx = get_exon_num(info.exons, info.cds.back()[1]);
    uint max_idx = get_exon_num(info.exons, cds_range["max"]);
    if (max_idx == last_idx) {
        info.cds.back()[1] = cds_range["max"];
    } else {
        info.cds.back()[1] = info.exons[last_idx][1];  // extend existing CDS
        std::vector<int> extra_cds = {info.exons[max_idx][0], cds_range["max"]};
        info.cds.push_back(extra_cds);
    }
}

// collect all features for a transcript into a single object
//
// When we load lines from gencode GTF files, each line represents a single exon
// or CDS, and we need to combine these based on transcript ID
void load_transcripts(std::vector<NamedTx> & transcripts, GTF &gtf_file, bool coding=true) {
    std::set<std::string> permit = {"exon", "CDS", "UTR", "transcript", 
        "stop_codon", "start_codon"};
    std::map<std::string, int> cds_range = {{"max", 0}, {"min", 999999999}};
    std::string tx_id = "";
    std::string symbol = "";
    std::string current;
    bool is_principal = false;
    TxInfo info;

    GTFLine gtf;

    while (true) {
        try {
            gtf = gtf_file.next();
        } catch (const std::out_of_range& e) {
            break;
        }
        

        if (permit.count(gtf.feature) == 0) {
            continue;
        } else if (coding and gtf.transcript_type != "protein_coding") {
            continue;
        }

        current = gtf.tx_id;
        if (tx_id == "") {
            tx_id = current;
            symbol = gtf.symbol;
        }

        if (tx_id != current) {
            // adjust CDS for start and stop codon coords
            include_end_codons(cds_range, info);
            Tx tx = Tx(info.name, info.chrom, info.start, info.end, info.strand[0]);
            tx.set_exons(info.exons, info.cds);
            tx.set_cds(info.cds);
            transcripts.push_back({symbol, tx, is_principal});
            info = {};
            tx_id = current;
            cds_range["max"] = 0;
            cds_range["min"] = 999999999;
            symbol = gtf.symbol;
            is_principal = false;
        }

        if (info.name == "") {
            info.name = tx_id;
            info.chrom = gtf.chrom;
            info.strand = gtf.strand;
        }

        if (gtf.is_principal) {
            is_principal = true;
        }

        if (gtf.feature == "transcript") {
            info.start = gtf.start;
            info.end = gtf.end;
        } else if (gtf.feature == "CDS") {
            info.cds.push_back(std::vector<int> {gtf.start, gtf.end});
            cds_range["max"] = std::max(std::max(cds_range["max"], gtf.start), gtf.end);
            cds_range["min"] = std::min(std::min(cds_range["min"], gtf.start), gtf.end);
        } else if (gtf.feature == "exon") {
            info.exons.push_back(std::vector<int> {gtf.start, gtf.end});
        } else if ((gtf.feature == "stop_codon") or (gtf.feature == "start_codon")) {
            cds_range["max"] = std::max(std::max(cds_range["max"], gtf.start), gtf.end);
            cds_range["min"] = std::min(std::min(cds_range["min"], gtf.start), gtf.end);
        }
    }

    // also include the final transcript (if it transcript exists)
    if (info.name != "") {
        include_end_codons(cds_range, info);
        Tx tx = Tx(info.name, info.chrom, info.start, info.end, info.strand[0]);
        tx.set_exons(info.exons, info.cds);
        tx.set_cds(info.cds);
        transcripts.push_back({symbol, tx, is_principal});
    }
}

std::vector<NamedTx> open_gencode(std::string path, bool coding) {
    GTF gtf_file(path);
    std::vector<NamedTx> transcripts;
    transcripts.reserve(80000);
    load_transcripts(transcripts, gtf_file, coding);
    return transcripts;
}

bool CompFunc(const GenePoint &l, const GenePoint &r) {
    return l.pos < r.pos;
}

std::vector<std::string> _in_region(std::string chrom, int start, int end, 
        std::map<std::string, std::vector<GenePoint>> & starts, 
        std::map<std::string, std::vector<GenePoint>> & ends,
        int max_window=2500000) {
    
    if (chrom.size() < 3 || chrom.substr(0, 3) != "chr") {
        chrom.insert(0, "chr");
    }
    
    if (starts.count(chrom) == 0) {
        throw std::invalid_argument("unknown_chrom: " + chrom);
    }
    
    std::vector<GenePoint> & chrom_starts = starts[chrom];
    std::vector<GenePoint> & chrom_ends = ends[chrom];
    std::set<std::size_t> inside;
    std::vector<std::string> symbols;
    symbols.reserve(std::max((end - start) / 50000, 1));  // expect 1 gene / 50 kb 
    
    // find indices to genes with a start inside the region
    int left_idx, right_idx;
    int idx = std::lower_bound(chrom_starts.begin(), chrom_starts.end(), GenePoint {start, "A"}, CompFunc) - chrom_starts.begin();
    left_idx = idx - 1;
    while (idx < (int) chrom_starts.size()) {
        GenePoint & edge = chrom_starts[idx];
        if (edge.pos > end) {
            break;
        } else {
            inside.insert(std::hash<std::string>{}(edge.symbol));
            symbols.push_back(edge.symbol);
        }
        idx += 1;
    }

    // find indices to genes with a end inside the region
    idx = (std::upper_bound(chrom_ends.begin(), chrom_ends.end(), GenePoint {end, "A"}, CompFunc) - chrom_ends.begin());
    right_idx = idx;
    idx = std::min(idx - 1, (int) chrom_ends.size() - 1);
    while (idx >= 0) {
        GenePoint & edge = chrom_ends[idx];
        if (edge.pos < start) {
            break;
        } else if (inside.count(std::hash<std::string>{}(edge.symbol)) == 0) {
            symbols.push_back(edge.symbol);
        }
        idx -= 1;
    }

    if (abs(end - start) > max_window) {
        // if the window is too wide to permit a gene to span it, just return
        return symbols;
    }
    // for genes that encapsulate the region, first find genes that start upstream
    static std::set<std::size_t> starts_before;
    starts_before.clear();
    for (; left_idx>=0; left_idx--) {
        GenePoint & edge = chrom_starts[left_idx];
        if (abs(edge.pos - end) > max_window) { // halt if distant from the region
            break;
        }
        starts_before.insert(std::hash<std::string>{}(edge.symbol));
    }

    // find genes that end downstream of the gene
    int length = (int) chrom_ends.size();
    for (; right_idx<length; right_idx++) {
        GenePoint & edge = chrom_ends[right_idx];
        if (abs(edge.pos - start) > max_window) { // halt if distant from the region
            break;
        }
        if (starts_before.count(std::hash<std::string>{}(edge.symbol)) != 0) {
            symbols.push_back(edge.symbol);
        }
    }

    return symbols;
}

} // namespace

// int main() {
//     std::string path = "/illumina/scratch/deep_learning/public_data/refdata/hg38/genes/gencode.v24.annotation.gtf";
//     gencode::open_gencode(path);
// }
// 
// g++ -std=c++11 gencode.cpp gtf.cpp tx.cpp gzstream/gzstream.C -Igzstream -lz



