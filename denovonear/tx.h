#ifndef DENOVONEAR_TX_H_
#define DENOVONEAR_TX_H_

#include <string>
#include <vector>

struct Region {
    int start;
    int end;
};

class Tx {
    std::string name;
    std::string chrom;
    int tx_start;
    int tx_end;
    char tx_strand;
    int cds_min;
    int cds_max;
    std::vector<Region> exons;
    std::vector<Region> cds;

 public:
    Tx(std::string transcript_id, std::string chromosome, int start_pos,
        int end_pos, char strand);
    void set_exons(std::vector<std::vector<int>> exon_ranges,
        std::vector<std::vector<int>> cds_ranges);
    void set_cds(std::vector<std::vector<int>> cds_ranges);
    Region fix_cds_boundary(int position);
    
    std::vector<Region> get_exons() { return exons; }
    std::vector<Region> get_cds() { return cds; }
    std::string get_name() { return name; }
    std::string get_chrom() { return chrom; }
    int get_start() { return tx_start; }
    int get_end() { return tx_end; }
    char get_strand() { return tx_strand; }
    int get_cds_start();
    int get_cds_end();
    
    bool in_exons(int position);
    Region find_closest_exon(int position);
    Region find_closest_exon(int position, std::vector<Region> ranges);
    bool in_coding_region(int position);
    int get_exon_containing_position(int position, std::vector<Region> ranges);
    int get_coding_distance(int pos_1, int pos_2);
    int chrom_pos_to_cds(int pos_1);
    
};

#endif  // DENOVONEAR_TX_H_
