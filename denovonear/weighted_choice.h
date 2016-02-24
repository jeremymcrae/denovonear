#ifndef DENOVONEAR_WEIGHTED_CHOICE_H_
#define DENOVONEAR_WEIGHTED_CHOICE_H_

#include <random>
#include <vector>

struct AlleleChoice {
    int pos;
    char ref;
    char alt;
};

class Chooser {
    std::vector<AlleleChoice> sites;
    std::vector<double> cumulative;
    std::uniform_real_distribution<double> dist;
    std::mt19937_64 generator;

 public:
    Chooser();
    void add_choice(int site, double prob, char ref='N', char alt='N');
    AlleleChoice choice();
    double get_summed_rate();
};

#endif  // DENOVONEAR_WEIGHTED_CHOICE_H_
