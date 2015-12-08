#include <random>
#include <vector>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <python2.7/Python.h>

#include "weighted_choice.h"
// g++ -std=c++0x -c -fPIC weighted_choice.cpp -o weighted_choice.o
// g++ -shared -Wl,-soname,libweightedchoice.so -o libweightedchoice.so weighted_choice.o

namespace weights {
    WeightedChoice::WeightedChoice()
    {
        /**
            Constructor for WeightedChoice class
        */
        
        // start the random sampler
        long long random_seed = std::chrono::system_clock::now().time_since_epoch().count();
        generator.seed(random_seed);
    }

    void WeightedChoice::add_choice(int site, double prob)
    {
         /**
            adds another choice to the class object
            
            @site site position (e.g. 100001)
            @prob site mutation rate (e.g. 0.000000005)
        */
        
        sites.push_back(site);
        
        double cumulative_sum = get_summed_rate();
        cumulative_sum += prob;
        cumulative.push_back(cumulative_sum);
        
        // any time we add another choice, reset the sampler, so we can sample
        // from all the possible entries.
        std::uniform_real_distribution<double> dist(0.0, cumulative_sum);
    }

    int WeightedChoice::choice()
    {
        /**
            chooses a random element using a set of probability weights
            
            @return the name of the randomly selected element (e.g. position)
        */
        
        // get a random float between 0 and the cumulative sum
        double number = dist(generator);
        
        // figure out where in the list a random probability would fall
        std::vector<double>::iterator pos;
        pos = std::lower_bound(cumulative.begin(), cumulative.end(), number);
        
        // return the site position matching the probability
        return sites[pos - cumulative.begin()];
    }

    double WeightedChoice::get_summed_rate()
    {
        /**
           gets the cumulative sum for all the current choices.
       */
       
        if (sites.empty()) {
            return 0.0;
        } else {
            return cumulative.back();
        }
    }
}
