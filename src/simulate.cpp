// #include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

#include "weighted_choice.h"

// gets the distances between all the pairs of elements from a list
//
// @param sites array of positions
// @return vector of paired positions
void _get_distances(int sites[], int & len, int distances[]) {
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

// check if any value in an vector is zero
//
// @param distances vector of values
// @return true/false for containing zero
bool _has_zero(int distances[], int & len) {
    return std::find(distances, distances + len, 0) != distances + len;
}

bool _has_zero(std::vector<int> distances) {
    int * dist = distances.data();
    int len = distances.size();
    return _has_zero(dist, len);
}

// gets the geometric mean of a vector of distances
//
// @param sites vector of distances
// @return geometric mean
double _geomean(int distances[], int & len) {
    bool zero_val = _has_zero(distances, len);
    double total = 0;
    // if some values are zero, adjust the values upwards, then add the log10
    // value, otherwise add the uncorrected log10 value
    if (zero_val) {
        for (int i=0; i < len; i++) {
            total += log(distances[i] + 1);
        }
    } else {
        for (int i=0; i < len; i++) {
            total += log(distances[i]);
        }
    }

    double mean = exp(total / len);
    // adjust mean back to where it should be if we had a zero value
    if (zero_val) { mean -= 1; }
    return mean;
}

double _geomean(std::vector<int> distances) {
    int * dist = distances.data();
    int len = distances.size();
    return _geomean(dist, len);
}

// simulates de novos weighted by mutation rate
//
// @param choices Chooser object, to sample sites
// @param iteration number of iterations to run
// @param de_novo_count number of de novos to simulate per iteration
// @return a list of mean distances for each iteration
std::vector<double> _simulate_distribution(Chooser & choices, int iterations,
    int de_novo_count) {
    
    // use a vector to return the mean distances, easier to call from python
    std::vector<double> mean_distances(iterations);
    
    int distance_len = (de_novo_count - 1) * (de_novo_count) / 2;
    int distances[distance_len];
    int positions[de_novo_count];
    // run through the required iterations
    for (int n=0; n < iterations; n++) {
        // randomly select de novo sites for the iteration
        for (int i=0; i < de_novo_count; i++) {
            positions[i] = choices.choice_pos_only();
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

// halt permutations if the P value could never be significant.
//
// assess whether the P-value could never fall below 0.1, and cut out after a
// smaller number of iterations, in order to minimise run time. Figure out the
// lower bound of the confidence interval for the current simulated P value.
//
// @param p_val current simulated P value
// @param iterations iterations run in order to obtain the simulated P value
// @param z standard normal deviate (eg 1.96 for 95% CI)
// @param alpha value above which to cease permuting
// @return True/False for whether to halt the permuations
bool _halt_permutation(double p_val, int iterations, double z, double alpha) {
    // @todo figure out whether this is legit, perhaps a Sequential
    // @todo Probability Ratio Test (SPRT) would be more appropriate.
    double delta = (z * sqrt((p_val * (1 - p_val))))/iterations;
    double lower_bound = p_val - delta;
    
    // if the lower bound of the confidence interval exceeds 0.1, then we
    // can be sure it's not going to ever get lower than 0.05.
    return lower_bound > alpha;
}

// simulates de novos weighted by mutation rate
//
// @param choices Chooser object, to sample sites
// @param iteration number of iterations to run
// @param de_novo_count number of de novos to simulate per iteration
// @param observed_value mean distance observed in the real de novo events
// @return a list of mean distances for each iteration
double _analyse_de_novos(Chooser & choices, int iterations, int de_novo_count,
    double observed_value) {
    
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
