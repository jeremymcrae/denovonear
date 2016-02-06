#ifndef DENOVONEAR_WEIGHTED_CHOICE_H_
#define DENOVONEAR_WEIGHTED_CHOICE_H_

#include <random>
#include <vector>

#define NULL nullptr

struct AlleleChoice {
    int pos;
    char ref;
    char alt;
}

class Chooser {
    std::vector<int> sites;
    std::vector<double> cumulative;
    std::vector< std::vector<char> > alleles;
    std::uniform_real_distribution<double> dist;
    std::mt19937_64 generator;

 public:
    Chooser();
    void add_choice(int site, double prob, char ref=NULL, char alt=NULL);
    int choice();
    AlleleChoice choice_with_alleles();
    double get_summed_rate();
};

#endif  // DENOVONEAR_WEIGHTED_CHOICE_H_
