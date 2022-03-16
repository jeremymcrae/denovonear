#ifndef DENOVONEAR_SITESCHECKS_H_
#define DENOVONEAR_SITESCHECKS_H_

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
    std::unordered_map<std::string, Chooser> rates;
    int boundary_dist;
    int kmer_length;
    int mid_pos;
    
    std::unordered_map<char, char> transdict = {
        {'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'}};
    
    std::vector<char> bases = {'A', 'C', 'G', 'T'};
    std::vector<std::string> categories = {"missense", "nonsense", "synonymous",
        "splice_lof", "splice_region", "loss_of_function", "intronic"};

 public:
    SitesChecks(Tx tx, std::vector<std::vector<std::string>> mut, bool cds_coords) :
         _tx { tx }, use_cds_coords { cds_coords } { init(mut); };
    SitesChecks(Tx tx, std::vector<std::vector<std::string>> mut, bool cds_coords, Tx mask) :
         _tx { tx }, masked { mask }, use_cds_coords { cds_coords } { has_mask = true; init(mut); };
    Chooser * __getitem__(std::string category) { return &rates[category]; };
    void initialise_choices();
    
    Tx _tx;
    void check_position(int bp);
    int get_offset(int bp);
    void check_consequence(std::string & cq, char & initial_aa, char & mutated_aa, int & offset);
    
 private:
    Tx masked = Tx("zz", "z", -100, -100, '+', "protein_coding");
    void init(std::vector<std::vector<std::string>> mut);
    bool has_mask = false;
    bool use_cds_coords = true;
};

Region _get_gene_range(Tx & tx);
char _get_mutated_aa(Tx & tx, char base, std::string codon, int intra_codon);

#endif  // DENOVONEAR_SITESCHECKS_H_
