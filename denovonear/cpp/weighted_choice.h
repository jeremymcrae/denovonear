#ifndef WEIGHTED_CHOICE_H
#define WEIGHTED_CHOICE_H

#include <random>
#include <vector>

namespace weights {
    class WeightedChoice
    {
        std::vector<int> sites;
        std::vector<double> cumulative;
        std::uniform_real_distribution<double> dist;
        std::mt19937_64 generator;

        public:
            WeightedChoice();
            void add_choice(int site, double prob);
            int choice();
            double get_summed_rate();
    };
}

#endif
