#ifndef DENOVONEAR_TX_H_
#define DENOVONEAR_TX_H_

#include <string>
#include <vector>
#include <unordered_map>

struct CDS_coords {
    int position;
    int offset;
};

struct Region {
    int start;
    int end;
};

struct Codon {
    int cds_pos;
    std::string codon_seq;
    int intra_codon;
    int codon_number;
    std::string initial_aa;
    int offset;
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
    std::string cds_sequence = "";
    int gdna_offset;
    std::string genomic_sequence = "";
    
    std::unordered_map<char, char> transdict = {
        {'a', 't'}, {'c', 'g'}, {'g', 'c'}, {'t', 'a'}, {'u', 'a'},
        {'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'}, {'U', 'A'}};
    
    std::unordered_map<std::string, std::string> aa_code = {
        {"AAA", "K"}, {"AAC", "N"}, {"AAG", "K"}, {"AAT", "N"},
        {"ACA", "T"}, {"ACC", "T"}, {"ACG", "T"}, {"ACT", "T"},
        {"AGA", "R"}, {"AGC", "S"}, {"AGG", "R"}, {"AGT", "S"},
        {"ATA", "I"}, {"ATC", "I"}, {"ATG", "M"}, {"ATT", "I"},
        {"CAA", "Q"}, {"CAC", "H"}, {"CAG", "Q"}, {"CAT", "H"},
        {"CCA", "P"}, {"CCC", "P"}, {"CCG", "P"}, {"CCT", "P"},
        {"CGA", "R"}, {"CGC", "R"}, {"CGG", "R"}, {"CGT", "R"},
        {"CTA", "L"}, {"CTC", "L"}, {"CTG", "L"}, {"CTT", "L"},
        {"GAA", "E"}, {"GAC", "D"}, {"GAG", "E"}, {"GAT", "D"},
        {"GCA", "A"}, {"GCC", "A"}, {"GCG", "A"}, {"GCT", "A"},
        {"GGA", "G"}, {"GGC", "G"}, {"GGG", "G"}, {"GGT", "G"},
        {"GTA", "V"}, {"GTC", "V"}, {"GTG", "V"}, {"GTT", "V"},
        {"TAA", "*"}, {"TAC", "Y"}, {"TAG", "*"}, {"TAT", "Y"},
        {"TCA", "S"}, {"TCC", "S"}, {"TCG", "S"}, {"TCT", "S"},
        {"TGA", "*"}, {"TGC", "C"}, {"TGG", "W"}, {"TGT", "C"},
        {"TTA", "L"}, {"TTC", "F"}, {"TTG", "L"}, {"TTT", "F"}};
    
    std::unordered_map<int, int> exon_to_cds;
    // std::map<int, int> exon_to_cds;
    void _cache_exon_cds_positions();
    
    void _fix_cds_length();

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
    Region find_closest_exon(int position, std::vector<Region> & ranges);
    bool in_coding_region(int position);
    int get_exon_containing_position(int position, std::vector<Region> & ranges);
    int get_coding_distance(int pos_1, int pos_2);
    CDS_coords chrom_pos_to_cds(int pos_1);
    
    int get_position_on_chrom(int cds_position, int offset=0);
    int get_codon_number_for_cds_position(int cds_position);
    int get_position_within_codon(int cds_position);
    void add_cds_sequence(std::string cds_dna);
    void add_genomic_sequence(std::string gdna, int offset);
    std::string get_cds_sequence() { return cds_sequence; }
    std::string get_genomic_sequence() { return genomic_sequence; }
    int get_genomic_offset() { return gdna_offset; }
    
    void _fix_transcript_off_by_one_bp();
    
    std::string reverse_complement(std::string seq);
    std::string get_centered_sequence(int pos, int length=3);
    std::string get_codon_sequence(int codon_number);
    std::string translate(std::string seq);
    
    Codon get_codon_info(int bp);
    int get_boundary_distance(int bp);
    
};

#endif  // DENOVONEAR_TX_H_
