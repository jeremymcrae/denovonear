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
    
    std::map<std::string, std::map<std::string, double>> mut_dict;
    std::map<std::string, Chooser> rates;
    int boundary_dist;
    int kmer_length;
    int mid_pos;
    
    std::map<std::string, std::string> transdict = {
        {"A", "T"}, {"T", "A"}, {"G", "C"}, {"C", "G"}};
    
    std::vector<std::string> bases = {"A", "C", "G", "T"};
    std::vector<std::string> categories = {"missense", "nonsense", "synonymous",
        "splice_lof", "splice_region", "loss_of_function"};

 public:
    SitesChecks(Tx tx, std::vector<std::vector<std::string>> mut) : _tx(tx) { init(mut); };
    SitesChecks(Tx tx, std::vector<std::vector<std::string>> mut, Tx mask) :
        _tx(tx), masked(mask) { has_mask = true; init(mut); };
    Chooser * __getitem__(std::string category) { return &rates[category]; };
    void initialise_choices();
    
    void check_position(int bp);
    std::string check_consequence(std::string initial_aa, std::string mutated_aa, int position);
    
 private:
    Tx _tx;
    Tx masked = Tx("zz", "z", -100, -100, '+');
    void init(std::vector<std::vector<std::string>> mut);
    bool has_mask = false;
};

Region _get_gene_range(Tx & tx);
std::string _get_mutated_aa(Tx & tx, std::string base, std::string codon, int intra_codon);

#endif  // DENOVONEAR_SITESCHECKS_H_
