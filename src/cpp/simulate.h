#ifndef SIMULATE_H
#define SIMULATE_H

std::vector<double> get_distances(int sites[], short len);
bool has_zero(std::vector<double> distances);
double get_geomean(std::vector<double> distances);
std::vector<double> simulate_distribution(WeightedChoice weights,
    int iterations, int de_novo_count);
double analyse_de_novos(int sites[], double probs[], int len,
    int iterations, int de_novo_count, double observed_value)

#endif