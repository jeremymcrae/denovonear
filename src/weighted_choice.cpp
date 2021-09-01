#include <random>
#include <vector>
#include <chrono>
#include <algorithm>

#include "weighted_choice.h"

Chooser::Chooser() {
    /**
        Constructor for Chooser class
    */
    
    // start the random sampler
    std::random_device rd;
    generator.seed(rd());
}

void Chooser::reset_sampler() {
    // any time we add another choice, reset the sampler, so we can sample
    // from all the possible entries.
    std::uniform_real_distribution<double> temp(0.0, get_summed_rate());
    dist = temp;
}

void Chooser::add_choice(int site, double prob, std::string ref, std::string alt, int offset) {
     /**
        adds another choice to the class object
        
        @site site position (e.g. 100001)
        @prob site mutation rate (e.g. 0.000000005)
        @ref reference allele for site e.g. 'A'
        @alt alternate allele for site e.g. 'T'
        @offset number of bases the position is offset from the true chromosomal
            site. This is relevant for non-CDS sites, such as within splice
            regions.
    */
    
    // keep track of the cumulative sum for each added site
    double cumulative_sum = get_summed_rate() + prob;
    cumulative.push_back(cumulative_sum);
    
    sites.push_back(AlleleChoice {site, ref, alt, prob, offset});
    reset_sampler();
}

// get index to randomly sampled element, using weighted probabilities
int Chooser::sampled_index() {
    if (cumulative.empty()) {
        return -1;
        // return AlleleChoice {-1, "N", "N", 0.0, 0};
    }
    
    // get a random float between 0 and the cumulative sum
    double number = dist(generator);
    
    // figure out where in the list a random probability would fall
    auto pos = std::lower_bound(cumulative.begin(), cumulative.end(), number);
    return pos - cumulative.begin();
}

// randomly sample a site, but only return the site position
// 
// This is a fast version for when we only want the site position, and want to
// avoid constructing a full AlleleChoice while returning.
int Chooser::choice_pos_only() {
    int idx = sampled_index();
    if (idx >= 0) {
        return sites[idx].pos;
    }
    return -1;
}

// randomly sample a site
AlleleChoice Chooser::choice() {
    int idx = sampled_index();
    if (idx >= 0) {
        return sites[idx];
    }
    return AlleleChoice {-1, "N", "N", 0.0, 0};
}

double Chooser::get_summed_rate() {
    /**
        gets the cumulative sum for all the current choices.
    */
    
    return (sites.empty()) ? 0.0 : cumulative.back() ;
}

void Chooser::append(Chooser other) {
    
    double current = get_summed_rate();
    int len = other.sites.size();
    for (int i=0; i < len; i++) {
        cumulative.push_back(other.cumulative[i] + current);
        sites.push_back(other.sites[i]);
    }
    
    reset_sampler();
}
