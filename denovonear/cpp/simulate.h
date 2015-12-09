#ifndef SIMULATE_H
#define SIMULATE_H

std::vector<double> get_distances(std::vector<int> sites);
bool has_zero(std::vector<double> distances);
double get_geomean(std::vector<double> distances);
bool halt_permutation(double p_val, int iterations, double z=10.0, double alpha=0.01);
std::vector<double> simulate_distribution(Chooser weights,
    int iterations, int de_novo_count);
double analyse_de_novos(Chooser weights, int iterations,
    int de_novo_count, double observed_value);

#endif
