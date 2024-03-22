#ifndef DENOVONEAR_SITESCHECKS_H_
#define DENOVONEAR_SITESCHECKS_H_

#include <cstdint>
#include <string>
#include <vector>
#include <map>

#include "tx.h"
#include "weighted_choice.h"

class SitesChecks {
    /**
    class to build weighted choice random samplers for nonsense, missense,
    and functional classes of variants, using site specific mutation rates
    
    Only include the site specific probability if it mutates to a different
    amino acid, or occurs close to an intron/exon boundary, using consequences
    defined at: http://www.ensembl.org/info/genome/variation/predicted_data.html
    */
    
    std::unordered_map<std::string, std::unordered_map<std::string, double>> mut_dict;
    std::unordered_map<std::uint64_t, std::unordered_map<char, double>> per_pos_rates;
    std::unordered_map<std::string, Chooser> rates;
    int boundary_dist;
    int kmer_length;
    int mid_pos;
    bool context_based=true;  // flag for whether the rates are based on per site 
                              // context, or from per-genome sites
    
    std::unordered_map<char, char> transdict = {
        {'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'}};
    
    std::vector<char> bases = {'A', 'C', 'G', 'T'};
    std::vector<std::string> categories = {"missense", "nonsense", "synonymous",
        "splice_lof", "splice_region", "loss_of_function", "intronic"};

 public:
    // this version uses sequence-context based rates
    SitesChecks(gencode::Tx tx, std::vector<std::vector<std::string>> mut, bool cds_coords) :
         _tx { tx }, use_cds_coords { cds_coords } { init(mut); };
    // this version uses rates fixed to each genome position
    SitesChecks(gencode::Tx tx, std::unordered_map<std::uint64_t, std::unordered_map<char, double>> mut, bool cds_coords) :
          per_pos_rates { mut }, _tx { tx }, use_cds_coords { cds_coords } { init(); };
    Chooser * __getitem__(std::string category) { return &rates[category]; };
    void initialise_choices();

    gencode::Tx _tx;
    void add_mask(gencode::Tx mask);
    void check_position(int bp);
    int get_offset(int bp);
    void check_consequence(std::string & cq, char & initial_aa, char & mutated_aa, int & offset);
    
 private:
    gencode::Tx masked = gencode::Tx("zz", "z", -100, -100, '+', "protein_coding");
    void init(std::vector<std::vector<std::string>> mut);
    void init();
    bool has_mask = false;
    bool use_cds_coords = true;
};

gencode::Region _get_gene_range(gencode::Tx & tx);
char _get_mutated_aa(gencode::Tx & tx, char base, std::string codon, int intra_codon);

#endif  // DENOVONEAR_SITESCHECKS_H_
