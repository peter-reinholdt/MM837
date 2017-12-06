#include <vector>
#include <time.h>
#include "toml.h"
#include <random>
#include <iostream>
#include "properties.h"

inline int delta_function(int i, int j);
inline void metropolis_sweep(std::vector<std::vector<int> >& lattice,
                             double beta,
                             int qmax,
                             std::mt19937& gen,
                             std::uniform_int_distribution<int>& q_dist, 
                             std::uniform_real_distribution<double>& double_dist);

inline void cluster_sweep(std::vector<std::vector<int> >& lattice,
                          double beta,
                          std::mt19937& gen,
                          std::uniform_int_distribution<int>& q_dist, 
                          std::uniform_int_distribution<int>& L_dist, 
                          std::uniform_real_distribution<double>& double_dist);


void metropolis(int n_steps_therm, int n_steps_prod, int side_length, int potts_q, double beta, std::string outfile, int outfreq, int conf_outfreq);
void cluster(int n_steps_therm, int n_steps_prod, int side_length, int potts_q, double beta, std::string outfile, int outfreq, int conf_outfreq);
void hybrid(int n_steps_therm, int n_steps_prod, int side_length, int potts_q, double beta, std::string outfile, int outfreq, int conf_outfreq, int metropolis_freq);

