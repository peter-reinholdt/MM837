#include <vector>
#include "toml.h"


inline void metropolis_sweep(){};
inline void compute_properties(){};


void metropolis(int n_steps_therm, int n_steps_prod, int side_length, int potts_q, double beta, std::string outfile, int outfreq){
    //store lattice as linear array, ordered as (i,j) --> side_length*i + j
    //65536 possible q values is probly enough.
    std::vector<unsigned short int> lattice;
    lattice.resize(side_length*side_length);
    
    //thermalize
    for (int n=0; n<n_steps_therm; n++){
        metropolis_sweep(lattice, side_length, potts_q, beta);
    }
    //production run
    for (int n=0; n<n_steps_prod; n++){
        metropolis_sweep();
        compute_properties();
    }
};



void cluster(int n_steps_therm, int n_steps_prod, int side_length, int potts_q, double beta, std::string outfile, int outfreq){};


