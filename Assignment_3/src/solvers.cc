#include <vector>
#include <time.h>
#include "toml.h"
#include <random>

inline void metropolis_sweep(std::vector<std::vector<int> >& lattice,
                             double beta,
                             std::uniform_int_distribution<int>& int_dist, 
                             std::uniform_real_distribution<float>& float_dist){
    int delta_E;
    int L = lattice.size();
    //first sweep through central part
    for (int i=1; i<L-1; i++){
        for (int j=1; j<L-1; j++){
            //fliiip
        }
    }
    //then sweep through borders
    //j=0 and j=side_length-1
    for (int i=0; i<L-1; i++){

    }
    //i=0 and i=side_length-1
    //only flip corners one time per sweep
    for (int j=1; j<L-2; j++){
    
    }
};
inline void compute_properties(){};


void metropolis(int n_steps_therm, int n_steps_prod, int side_length, int potts_q, double beta, std::string outfile, int outfreq){
    //store lattice as 2d array
    //65536 possible q values is probly enough.
    std::vector<std::vector<int> > lattice;
    std::mt19937 gen(time(NULL));
    std::uniform_int_distribution<int> int_dist(0, potts_q-1);
    std::uniform_real_distribution<float> real_dist(0,1);
    
    //thermalize
    for (int n=0; n<n_steps_therm; n++){
        metropolis_sweep(lattice, beta, int_dist, real_dist);
    }
    //production run
    for (int n=0; n<n_steps_prod; n++){
        //metropolis_sweep();
        //compute_properties();
    }
};



void cluster(int n_steps_therm, int n_steps_prod, int side_length, int potts_q, double beta, std::string outfile, int outfreq){};


