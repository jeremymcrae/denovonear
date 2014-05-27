#include <iostream>
#include <random>
#include <string>
#include <array>
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
    std::mt19937_64 generator;
    std::vector<double> cumulative;
    double cumulative_sum;
    
    public:
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
    int random_seed;
    random_seed = std::chrono::system_clock::now().time_since_epoch().count();
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

std::vector<double> get_distances(int sites[], short len)
{
    /**
        makes a vector of all the two element combinations of two elements
        
        @sites array of positions
        @len length of the sites array
        @return vector of paired positions
    */
    
    std::vector<double> distances;
    
    // get all non-repeating combinations of the sites
    for (int i=0; i<len; i++)
    {
        for (int j=i; j<len; j++)
        {
            // only include if the array positions differ, so we avoid finding
            // the distance to itself
            if (j != i)
            {
                double distance = abs(sites[i] - sites[j]);
                distances.push_back(distance);
            }            
        }
    }
    
    return distances;
}

bool has_zero(std::vector<double> distances)
{
    /**
        @check if any value in a vector is zero
        
        @distances vector of values
        @return true/false for containing zero
    */
        
    unsigned sz = distances.size();
    bool zero_val = false;
    
    // find if any of the elements are zero
    for (unsigned i=0; i<sz; i++)
    {
        if (distances[i] == 0.0) 
        {
            zero_val = true;
            break;
        }
    }
    
    return zero_val;
}

double get_geomean(std::vector<double> distances)
{
    /**
        gets the geometric mean of a vector of distances
        
        @sites array of positions
        @return geometric mean
    */
    unsigned sz = distances.size();
    bool zero_val = has_zero(distances);
    
    if (zero_val) 
    {
        // if some values are zero, adjust all of the values upwards, and get 
        // the log10 value
        for (unsigned i=0; i<sz; i++) distances[i] = log10(distances[i] + 1);
    }
    else 
    {
        // get the log10 value when we lack zero values
        for (unsigned i=0; i<sz; i++) distances[i] = log10(distances[i]);
    }
    
    // sum the distances in the vector
    double total = 0;
    for (unsigned i=0; i<sz; i++) total += distances[i];
    
    // calculate the mean value
    double mean = total/sz;
    mean = std::pow(10, mean);
    
    // adjust mean back to where it should be if we had a zero value
    if (zero_val) mean -= 1;
    
    return mean;
}

PyObject* simulate_distribution(int sites[], double probs[], int len,
    int iterations, int de_novo_count)
{
    /**
        simulates de novos weighted by mutation rate
        
        @sites array of CDS positions
        @probs array of mutation rates matches to the sites
        @len length of the sites and probs arrays
        @iteration number of iterations to run
        @de_novo_count number of de novos to simulate per iteration
        @return a list of mean distances for each iteration
    */
    
    // use a python object to return the mean distances, makes it easier to
    // call from python
    PyObject* mean_distances = PyList_New(0);
    
    // construct the weighted sampler
    WeightedChoice weights (sites, probs, len);
    
    // run through the required iterations
    for (int n=0; n < iterations; n++)
    {
        // randomly select de novo sites for the iteration
        int positions[de_novo_count];
        for (int i=0; i < de_novo_count; i++) 
        {
            int choice = weights.choice();
            positions[i] = choice;
        }
        
        // convert the positions into distances between all pairs, and get the
        // geometric mean distance of all the distances
        std::vector<double> distances = get_distances(positions, de_novo_count);
        double mean_distance = get_geomean(distances);
        PyList_Append(mean_distances, PyFloat_FromDouble(mean_distance));
    }
    
    return mean_distances;
}

extern "C" {
    WeightedChoice* WeightedChoice_new(int pos[], double p[], int len) { return new WeightedChoice(pos, p, len); }
    int WeightedChoice_choice(WeightedChoice* chooser) { return  chooser->choice(); }
    PyObject* c_simulate_distribution(int pos[], double p[], int len, 
        int iterations, int de_novo_count) 
    { 
        return simulate_distribution(pos, p, len, iterations, de_novo_count); 
    }
}
