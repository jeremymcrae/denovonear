#ifndef CHOOSER_H
#define CHOOSER_H

#include <random>
#include <vector>

class Chooser
{
    std::vector<int> sites;
    std::vector<double> cumulative;
    std::uniform_real_distribution<double> dist;
    std::mt19937_64 generator;

    public:
        Chooser();
        void add_choice(int site, double prob);
        int choice();
        double get_summed_rate();
};

#endif
