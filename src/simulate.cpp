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
    int * sites_array = new int[size];
    for (int i=0; i<size; i++) {
        sites_array[i] = sites[i];
    }
    int * dist_array = new int[len];
    _get_distances(sites_array, size, dist_array);
    distances.resize(len);
    
    for (int i=0; i<len; i++) {
        distances[i] = dist_array[i];
    }
    delete[] sites_array;
    delete[] dist_array;
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
double _geomean_small(int distances[], int & len) {
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

// gets the geometric mean of a vector of distances
//
// This is optimized for large arrays, it is faster than the log-transform
// based method if there are >300 items in the array. This occurs with >17
// de novo mutations for checking pairwise distances.
//
// @param distances array of distances
// @param len number of items in array
// @return geometric mean
double _geomean_large(int distances[], int & len) {
    // adapted from:
    // https://stackoverflow.com/questions/19980319/efficient-way-to-compute-geometric-mean-of-many-numbers
    bool zero_val = _has_zero(distances, len);
    if (zero_val) {
        for (int i=0; i<len; i++) {
            distances[i] += 1;
        }
    }
    long long ex = 0;
    auto do_bucket = [&distances, &ex](int first, int last) -> double
    {
        double ans = 1.0;
        for ( ;first != last; ++first) {
            int i;
            ans *= std::frexp(static_cast<double>(distances[first]), &i);
            ex += i;
        }
        return ans;
    };

    const int bucket_size = (int)-std::log2(std::numeric_limits<double>::min());
    std::size_t buckets = len / bucket_size;

    double invN = 1.0 / len;
    double m = 1.0;

    for (std::size_t i = 0;i < buckets; ++i) {
        m *= std::pow(do_bucket(i * bucket_size, (i+1) * bucket_size), invN);
    }

    m *= std::pow(do_bucket( buckets * bucket_size, len), invN);
    double mean = std::pow( std::numeric_limits<double>::radix, ex * invN ) * m;
    if (zero_val) { mean -= 1; }
    return mean;
}

// get the geometric mean of an array of values
//
double _geomean(int distances[], int & len) {
    if (len < 300) {
        return _geomean_small(distances, len);
    } else {
        return _geomean_large(distances, len);
    }
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
    int * distances = new int[distance_len];
    int * positions = new int[de_novo_count];
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
    delete[] distances;
    delete[] positions;
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

    std::uint32_t n_smaller = 0;
    while (iterations < 100000000 && sim_prob == minimum_prob) {
        int iters_to_run = iterations - dist.size();
        
        minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));
        
        // simulate mean distances between de novos
        dist = _simulate_distribution(choices, iters_to_run, de_novo_count);
        
        for (auto &x : dist) {
            n_smaller += (x <= observed_value);
        }
        
        // estimate the probability from the number of times a random distance
        // is equal to or smaller than the observed value
        sim_prob = (1.0 + (double)n_smaller) / (1.0 + (double)iterations);
        
        // halt permutations if the P value could never be significant
        double z = 10.0;
        double alpha = 0.1;
        if (_halt_permutation(sim_prob, iterations, z, alpha)) { break; }
        
        iterations += 1000000;  // for if we need to run more iterations
    }
    
    return sim_prob;
}
