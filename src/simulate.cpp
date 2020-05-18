// #include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

#include "weighted_choice.h"

void _get_distances(int sites[], int & len, int distances[]) {
    /**
        gets the distances between all the pairs of elements from a list
        
        @sites array of positions
        @return vector of paired positions
    */
    // get all non-repeating combinations of the sites
    int idx = 0;
    for (int i=0; i < len; i++) {
        for (int j=i+1; j < len; j++) {
            // only include if the array positions differ, so we avoid finding
            // the distance to itself
            distances[idx] = abs(sites[i] - sites[j]);
            idx += 1;
        }
    }
}

void _get_distances(std::vector<int> sites, std::vector<int> & distances) {
    int size = sites.size();
    int len = ((size - 1) * size) / 2;
    int sites_array[size];
    for (int i=0; i<size; i++) {
        sites_array[i] = sites[i];
    }
    int dist_array[len];
    _get_distances(sites_array, size, dist_array);
    distances.resize(len);
    
    for (int i=0; i<len; i++) {
        distances[i] = dist_array[i];
    }
}

bool _has_zero(int distances[], int & len) {
    /**
        @check if any value in an vector is zero
        
        @distances vector of values
        
        @return true/false for containing zero
    */
    return std::find(distances, distances + len, 0) != distances + len;
}

bool _has_zero(std::vector<int> distances) {
    int * dist = distances.data();
    int len = distances.size();
    return _has_zero(dist, len);
}

double _geomean(int distances[], int & len) {
    /**
        gets the geometric mean of a vector of distances
        
        @sites vector of distances
        @return geometric mean
    */
    bool zero_val = _has_zero(distances, len);
    
    double total = 0;
    // if some values are zero, adjust the values upwards, then add the log10
    // value, otherwise add the uncorrected log10 value
    if (zero_val) {
        for (int i=0; i < len; i++) {
            total += log10(distances[i] + 1);
        }
    } else {
        for (int i=0; i < len; i++) {
            total += log10(distances[i]);
        }
    }
    
    // calculate the mean value
    double mean = total/len;
    mean = std::pow(10, mean);
    
    // adjust mean back to where it should be if we had a zero value
    if (zero_val) { mean -= 1; }
    
    return mean;
}

double _geomean(std::vector<int> distances) {
    int * dist = distances.data();
    int len = distances.size();
    return _geomean(dist, len);
}

std::vector<double> _simulate_distribution(Chooser & choices, int iterations,
    int de_novo_count) {
    /**
        simulates de novos weighted by mutation rate
        
        @choices Chooser object, to sample sites
        @iteration number of iterations to run
        @de_novo_count number of de novos to simulate per iteration
        @return a list of mean distances for each iteration
    */
    
    // use a vector to return the mean distances, easier to call from python
    std::vector<double> mean_distances(iterations);
    
    int distance_len = (de_novo_count - 1) * (de_novo_count) / 2;
    int distances[distance_len];
    int positions[de_novo_count];
    // run through the required iterations
    for (int n=0; n < iterations; n++) {
        // randomly select de novo sites for the iteration
        for (int i=0; i < de_novo_count; i++) {
            positions[i] = choices.choice().pos;
        }
        
        // convert the positions into distances between all pairs, and get the
        // geometric mean distance of all the distances
        _get_distances(positions, de_novo_count, distances);
        
        mean_distances[n] = _geomean(distances, distance_len);
    }
    
    // make sure the mean distances are sorted, so we can quickly merge with
    // previous distances
    std::sort(mean_distances.begin(), mean_distances.end());
    
    return mean_distances;
}

bool _halt_permutation(double p_val, int iterations, double z, double alpha) {
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
    return lower_bound > alpha;
}

double _analyse_de_novos(Chooser & choices, int iterations, int de_novo_count,
    double observed_value) {
    /**
        simulates de novos weighted by mutation rate
        
        @choices Chooser object, to sample sites
        @iteration number of iterations to run
        @de_novo_count number of de novos to simulate per iteration
        @observed_value mean distance observed in the real de novo events
        @return a list of mean distances for each iteration
    */
    
    double minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));
    double sim_prob = minimum_prob;
    std::vector<double> dist;
    
    while (iterations < 100000000 && sim_prob == minimum_prob) {
        int iters_to_run = iterations - dist.size();
        
        minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));
        
        // simulate mean distances between de novos
        std::vector<double> new_dist = _simulate_distribution(choices,
            iters_to_run, de_novo_count);
        
        // merge the two sorted lists into a sorted vector
        std::vector<double> v(iterations);
        std::merge(dist.begin(), dist.end(), new_dist.begin(), new_dist.end(),
            v.begin());
        dist = v;
        
        // figure out where in the list a random probability would fall
        std::vector<double>::iterator pos;
        pos = std::upper_bound(dist.begin(), dist.end(), observed_value);
        double position = pos - dist.begin();
        
        // estimate the probability from the position
        sim_prob = (1.0 + position)/(1.0 + dist.size());
        
        // halt permutations if the P value could never be significant
        double z = 10.0;
        double alpha = 0.1;
        if (_halt_permutation(sim_prob, iterations, z, alpha)) { break; }
        
        iterations += 1000000;  // for if we need to run more iterations
    }
    
    return sim_prob;
}
