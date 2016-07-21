#ifndef DENOVONEAR_SITERATES_H_
#define DENOVONEAR_SITERATES_H_

#include <string>
#include <vector>
#include <map>

#include <tx.h>
#include <weighted_choice.h>

// Region get_gene_range(Tx tx);
// std::string get_mutated_aa(Tx tx, std::string base, std::string codon, int intra_codon);

class SitesChecks {
    /**
    class to build weighted choice random samplers for nonsense, missense,
    and functional classes of variants, using site specific mutation rates
    
    Only include the site specific probability if it mutates to a different
    amino acid, or occurs close to an intron/exon boundary, using consequences
    defined at: http://www.ensembl.org/info/genome/variation/predicted_data.html
    */
    
    Tx _tx;
    std::map<std::string std::map<std::string double>> mut_dict;
    Tx masked;
    int boundary_dist;
    
    std::map<char, char> transdict = {
        {'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'};
    
    std::vector<std::string> bases = {"A", "C", "T", "G"};
    std::vector<std::string> categories = {"missense", "nonsense", "synonymous",
        "splice_lof", "splice_region", "loss_of_function"};

 public:
    SitesChecks(Tx tx, std::vector<std::vector<std::string>> rates, Tx masked_sites);
    WeightedChoice __getitem__(std::string) { return rates[category]; };
    bool splice_lof_check(std::string initial_aa, std::string mutated_aa, int position);
    bool nonsense_check(std::string initial_aa, std::string mutated_aa, int position);
    bool missense_check(std::string initial_aa, std::string mutated_aa, int position);
    bool splice_region_check(std::string initial_aa, std::string mutated_aa, int position);
    bool synonymous_check(std::string initial_aa, std::string mutated_aa, int position);
    
    void check_position(int bp);
    
};

#endif  // DENOVONEAR_SITERATES_H_
