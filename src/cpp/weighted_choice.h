#ifndef WEIGHTED_CHOICE_H
#define WEIGHTED_CHOICE_H

#include <random>
#include <vector>

class WeightedChoice 
{
    int * sites;
    double * probs;
    int length;
    std::vector<double> cumulative;
    double cumulative_sum;
    
    public:
        std::mt19937_64 generator;
        WeightedChoice (int sites[], double probs[], int len);
        std::vector<double> make_cumulative_sums();
        int choice();
};

#endif