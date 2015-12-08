#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <python2.7/Python.h>
#include "weighted_choice.h"


double* get_distances(int sites[], short len, int combos)
{
    /**
        gets the distances between all the pairs of elements from a list
        
        @sites array of positions
        @len length of the sites array
        @combos number of pairwise combinations possible from the sites
        @return vector of paired positions
    */
    
    double* distances = new double [combos];
    
    int pos = 0;
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
                distances[pos] = distance;
                pos += 1;
            }
        }
    }
    
    return distances;
}

bool has_zero(double distances[], int sz)
{
    /**
        @check if any value in an array is zero
        
        @distances array of values
        @sz length of the distances array
        @return true/false for containing zero
    */
    
    bool zero_val = false;
    
    // find if any of the elements are zero
    for (int i=0; i<sz; i++)
    {
        if (distances[i] == 0.0)
        {
            zero_val = true;
            break;
        }
    }
    
    return zero_val;
}

double get_geomean(double distances[], int sz)
{
    /**
        gets the geometric mean of a vector of distances
        
        @sites array of positions
        @sz length of the distances array
        @return geometric mean
    */
    
    bool zero_val = has_zero(distances, sz);
    
    if (zero_val)
    {
        // if some values are zero, adjust all of the values upwards, and get
        // the log10 value
        for (int i=0; i<sz; i++) distances[i] = log10(distances[i] + 1);
    }
    else
    {
        // get the log10 value when we lack zero values
        for (int i=0; i<sz; i++) distances[i] = log10(distances[i]);
    }
    
    // sum the distances in the vector
    double total = 0;
    for (int i=0; i<sz; i++) total += distances[i];
    
    // calculate the mean value
    double mean = total/sz;
    mean = std::pow(10, mean);
    
    // adjust mean back to where it should be if we had a zero value
    if (zero_val) mean -= 1;
    
    return mean;
}

std::vector<double> simulate_distribution(WeightedChoice weights, long iterations,
    int de_novo_count)
{
    /**
        simulates de novos weighted by mutation rate
        
        @weights WeightedChoice object, to sample sites
        @iteration number of iterations to run
        @de_novo_count number of de novos to simulate per iteration
        @return a list of mean distances for each iteration
    */
    
    // construct the weighted sampler
    // TODO: figure out why I can't start this in the function that calls this,
    // TODO: when I tried, each set of simulations gave the same choices, which
    // TODO: implies the sampler restarts each time with the same seed.
    // WeightedChoice weights();
    
    // use a vector to return the mean distances, easier to call from python
    int combos = ((de_novo_count - 1) * de_novo_count)/2;
    std::vector<double> mean_distances;
    
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
        double *distances = get_distances(positions, de_novo_count, combos);
        double mean_distance = get_geomean(distances, combos);
        delete[] distances;
        mean_distances.push_back(mean_distance);
    }
    
    // make sure the mean distances are sorted, so we can quickly merge with
    // previous distances
    std::sort(mean_distances.begin(), mean_distances.end());
    
    return mean_distances;
}

bool halt_permutation(double p_val, int iterations, double z, double alpha)
{
    /**
        halt permutations if the P value could never be significant
    
        assess whether the P-value could never fall below 0.1, and cut
        out after a smaller number of iterations, in order to minimise
        run time. Figure out the lower bound of the confidence interval
        for the current simulated P value.
        TODO: figure out whether this is legit, perhaps a Sequential
        TODO: Probability Ratio Test (SPRT) would be more appropriate.
        
        @p_val current simulated P value
        @iterations iterations run in order to obtain the simulated P value
        @z standard normal deviate (eg 1.96 for 95% CI)
        @alpha value above which to cease permuting
        @return True/False for whether to halt the permuations
    */
    double delta = (z * sqrt((p_val * (1 - p_val))))/iterations;
    double lower_bound = p_val - delta;
    
    // if the lower bound of the confidence interval exceeds 0.1, then we
    // can be sure it's not going to ever get lower than 0.05.
    bool exceeds = lower_bound > alpha;
    
    return exceeds;
}

double analyse_de_novos(WeightedChoice weights, int iterations, int de_novo_count,
    double observed_value)
{
    /**
        simulates de novos weighted by mutation rate
        
        @weights WeightedChoice object, to sample sites
        @iteration number of iterations to run
        @de_novo_count number of de novos to simulate per iteration
        @observed_value mean distance observed in the real de novo events
        @return a list of mean distances for each iteration
    */
    
    double minimum_prob = 1.0/(1.0 + (double) iterations);
    double sim_prob = minimum_prob;
    std::vector<double> dist;
    
    while (iterations < 100000000 and sim_prob == minimum_prob)
    {
        unsigned long iters_to_run = iterations - dist.size();
        
        minimum_prob = 1.0/(1.0 + (double) iterations);
        
        // simulate mean distances between de novos
        std::vector<double> new_dist = simulate_distribution(weights, iters_to_run, de_novo_count);
        
        // merge the two sorted lists into a sorted vector
        std::vector<double> v(iterations);
        std::merge(dist.begin(), dist.end(), new_dist.begin(), new_dist.end(), v.begin());
        dist = v;
        
        // figure out where in the list a random probability would fall
        std::vector< double >::iterator pos;
        pos = std::upper_bound(dist.begin(), dist.end(), observed_value);
        double position = pos - dist.begin();
        
        // estimate the probability from the position
        sim_prob = (1.0 + position)/(1.0 + dist.size());
        
        // halt permutations if the P value could never be significant
        double z = 10.0;
        double alpha = 0.1;
        if (halt_permutation(sim_prob, iterations, z, alpha)) { break; }
        
        iterations += 100000; // for if we need to run more iterations
    }
    
    return sim_prob;
}


extern "C" {
    WeightedChoice* WeightedChoice_new() {
        return new WeightedChoice();
    }
    void WeightedChoice_add_choice(WeightedChoice* chooser, int site, double prob) {
        return chooser->add_choice(site, prob);
    }
    int WeightedChoice_choice(WeightedChoice* chooser) {
        return chooser->choice();
    }
    double WeightedChoice_get_summed_rate(WeightedChoice* chooser) {
        return chooser->get_summed_rate();
    }
}
