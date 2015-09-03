// #include <iostream>
#include <random>
// #include <string>
// #include <array>
#include <vector>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <python2.7/Python.h>


// g++ -std=c++0x -c -fPIC weighted_choice.cpp -o weighted_choice.o
// g++ -shared -Wl,-soname,libweightedchoice.so -o libweightedchoice.so weighted_choice.o

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

WeightedChoice::WeightedChoice(int pos[], double p[], int len)
{
    /**
        Constructor for WeightedChoice class
        
        @sites array of variant positions
        @p array of variant mutation rates
        @len length of the sites array
    */
    
    sites = pos;
    
    probs = p;
    length = len;
    
    // start the random sampler
    long long random_seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed(random_seed);
    cumulative = make_cumulative_sums();
}

std::vector< double > WeightedChoice::make_cumulative_sums()
{
     /**
        generates an array of cumulative probabilities
        
        @return array of cumulative probabilities corresponding to each element
            of the variants array
    */
    
    std::vector< double > cumulative;
    cumulative_sum = 0.0;
    for (int n=0; n<length; n++)
    {
        cumulative_sum += probs[n];
        cumulative.push_back(cumulative_sum);
    }
    
    return cumulative;
}

int WeightedChoice::choice()
{
    /**
        chooses a random element using a set of probability weights
        
        @return the name of the randomly selected element (e.g. position)
    */
    
    // get a random float between 0 and the cumulative sum
    std::uniform_real_distribution<double> dist(0.0, cumulative_sum);
    double number = dist(generator);
    
    // figure out where in the list a random probability would fall
    std::vector< double >::iterator pos;
    pos = std::lower_bound(cumulative.begin(), cumulative.end(), number);
    
    int position = sites[pos - cumulative.begin()];
    
    // return the site position matching the probability
    return position;
}
