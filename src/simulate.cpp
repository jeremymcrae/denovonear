// #include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <stdexcept>

#include "weighted_choice.h"
#include "simulate.h"

// gets the CDS distances between all the pairs of CDS positions
//
// @param sites array of positions
// @param distances array of pairwise distances to fill in
void _get_distances(std::vector<int> & sites, std::vector<int> & distances) {
    int size = sites.size();
    int len = ((size - 1) * size) / 2;
    distances.resize(len);
    // get all non-repeating combinations of the sites
    int idx = 0;
    for (int i=0; i < size; i++) {
        for (int j=i+1; j < size; j++) {
            // only include if the array positions differ, so we avoid finding
            // the distance to itself
            distances[idx] = abs(sites[i] - sites[j]);
            idx += 1;
        }
    }
}

// gets the 3D distances between all the pairs of amino acid coordinates
//
// @param sites vector of coords
// @param distances vector of pairwise distances to fill in
void _get_structure_distances(std::vector<Coord> & sites, std::vector<double> & distances) {
    int size = sites.size();
    int len = ((size - 1) * size) / 2;
    distances.resize(len);
    
     // get all non-repeating combinations of the sites
    double x_delta, y_delta, z_delta;
    int idx = 0;
    for (int i=0; i < size; i++) {
        for (int j=i+1; j < size; j++) {
            // only include if the array positions differ, so we avoid finding
            // the distance to itself
            x_delta = (sites[i].x - sites[j].x);
            y_delta = (sites[i].y - sites[j].y);
            z_delta = (sites[i].z - sites[j].z);
            
            x_delta = x_delta * x_delta;
            y_delta = y_delta * y_delta;
            z_delta = z_delta * z_delta;
            
            distances[idx] = std::sqrt(x_delta + y_delta + z_delta);
            idx += 1;
        }
    }
}

// check if any value in an vector is zero
//
// @param distances vector of values
// @return true/false for containing zero
template<typename T>
bool _has_zero(std::vector<T> & distances) {
    return std::find(distances.begin(), distances.end(), 0) != distances.end();
}

// gets the geometric mean of a vector of distances
//
// @param sites vector of distances
// @return geometric mean
template<typename T>
double _geomean_small(std::vector<T> & distances) {
    bool zero_val = _has_zero(distances);
    double total = 0;
    int len = distances.size();
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
// @param distances vector of distances
// @return geometric mean
template<typename T>
double _geomean_large(std::vector<T> & distances) {
    // adapted from:
    // https://stackoverflow.com/questions/19980319/efficient-way-to-compute-geometric-mean-of-many-numbers
    bool zero_val = _has_zero(distances);
    int len = distances.size();
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

// get the geometric mean of an vector of ints
double _geomean(std::vector<int> & distances) {
    if (distances.size() < 300) {
        return _geomean_small(distances);
    } else {
        return _geomean_large(distances);
    }
}

// get the geometric mean of an vector of doubles
double _geomean_double(std::vector<double> & distances) {
    if (distances.size() < 300) {
        return _geomean_small(distances);
    } else {
        return _geomean_large(distances);
    }
}

// simulates CDS distances between de novos weighted by mutation rate
//
// @param choices Chooser object, to sample sites
// @param iteration number of iterations to run
// @param de_novo_count number of de novos to simulate per iteration
// @return a list of mean distances for each iteration
std::vector<double> _simulate_distances(Chooser & choices, int iterations,
    int de_novo_count) {
    
    // use a vector to return the mean distances, easier to call from python
    std::vector<double> mean_distances(iterations);
    
    int distance_len = (de_novo_count - 1) * (de_novo_count) / 2;
    std::vector<int> distances(distance_len);
    std::vector<int> positions(de_novo_count);
    // run through the required iterations
    for (int n=0; n < iterations; n++) {
        // randomly select de novo sites for the iteration
        for (int i=0; i < de_novo_count; i++) {
            positions[i] = choices.choice_pos_only();
        }
        
        // convert the positions into distances between all pairs, and get the
        // geometric mean distance of all the distances
        _get_distances(positions, distances);
        
        mean_distances[n] = _geomean(distances);
    }
    return mean_distances;
}

// simulates distances between amino acids for de novo mutations
//
// @param choices Chooser object, to sample sites
// @param coords vector of 3D coords
// @param iteration number of iterations to run
// @param de_novo_count number of de novos to simulate per iteration
// @return a list of mean distances for each iteration
std::vector<double> _simulate_structure_distances(Chooser & choices, 
                                                  std::vector<Coord> & coords, 
                                                  int iterations,
                                                  int de_novo_count) {
    
    if (((int) coords.size() * 9) + 6 < choices.len()) {
        throw std::invalid_argument("number of coords is too low for the transcript site sampler");
    }
    
    // use a vector to return the mean distances, easier to call from python
    std::vector<double> mean_distances(iterations);
    
    int distance_len = (de_novo_count - 1) * (de_novo_count) / 2;
    std::vector<double> distances(distance_len);
    std::vector<Coord> positions(de_novo_count);
    // run through the required iterations
    for (int n=0; n < iterations; n++) {
        // randomly select de novo sites for the iteration
        
        std::size_t idx;
        for (int i=0; i < de_novo_count; i++) {
            idx = choices.choice_pos_only() / 3;
            while (idx >= coords.size()) {
                // ensure we don't use the terminator codon
                idx = choices.choice_pos_only() / 3;
            }
            positions[i] = coords[idx];
        }
        
        // convert the positions into distances between all pairs, and get the
        // geometric mean distance of all the distances
        _get_structure_distances(positions, distances);
        
        mean_distances[n] = _geomean_double(distances);
    }
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

// simulates clustering of de novos weighted by mutation rate
//
// @param choices Chooser object, to sample sites
// @param iteration number of iterations to run
// @param de_novo_count number of de novos to simulate per iteration
// @param observed_value mean distance observed in the real de novo events
// @return a list of mean distances for each iteration
double _simulate_clustering(Chooser & choices, int iterations, int de_novo_count,
    double observed_value) {
    
    double minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));
    double sim_prob = minimum_prob;
    std::vector<double> dist;

    std::uint32_t n_smaller = 0;
    while (iterations < 100000000 && sim_prob == minimum_prob) {
        int iters_to_run = iterations - dist.size();
        
        minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));
        
        // simulate mean distances between de novos
        dist = _simulate_distances(choices, iters_to_run, de_novo_count);
        
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

// simulates 3D clustering of de novos weighted by mutation rate
//
// @param choices Chooser object, to sample sites
// @param coords vector of 3D coords
// @param iteration number of iterations to run
// @param de_novo_count number of de novos to simulate per iteration
// @param observed_value mean distance observed in the real de novo events
// @return a list of mean distances for each iteration
double _simulate_structure_clustering(Chooser & choices, 
                                      std::vector<Coord> & coords,
                                      int iterations,
                                      int de_novo_count,
                                      double observed_value) {
    
    double minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));
    double sim_prob = minimum_prob;
    std::vector<double> dist;

    std::uint32_t n_smaller = 0;
    while (iterations < 100000000 && sim_prob == minimum_prob) {
        int iters_to_run = iterations - dist.size();
        
        minimum_prob = 1.0/(1.0 + static_cast<double>(iterations));
        
        // simulate mean distances between de novos
        dist = _simulate_structure_distances(choices, coords, iters_to_run, de_novo_count);
        
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
