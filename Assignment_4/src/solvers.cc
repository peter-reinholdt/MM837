#include <vector>
#include "toml.h"
#include <random>
#include <iostream>
#include <cmath>



inline void to_interval(double& x){
    //-PI..PI interval
    x = fmod(x, M_PI);
}

void lattice_to_interval(std::vector<std::vector<double> >& lattice){
    int L = lattice.size();
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            to_interval(lattice[i][j];
        }
    }
}


inline void metropolis_sweep(std::vector<std::vector<double> >& lattice, 
                             xoroshiro128plus& gen, 
                             std::uniform_real_distribution<double>& delta_dist){
    int L = lattice.size();
    int i_plus, j_plus, i_minus, j_minus;
    double delta_E;

    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            delta_E = 0.0;
            i_plus  = (i != l-1) ? i+1: 0;
            i_minus = (i != 0  ) ? i-1: l-1;
            j_plus  = (j != l-1) ? j+1: 0;
            j_minus = (j != 0  ) ? j-1: l-1;
            
            //neighboring angles
            n0 = lattice[i_plus]][j];
            n1 = lattice[i_minus]][j];
            n2 = lattice[i]][j_plus];
            n3 = lattice[i]][j_minus];

            //delta_E = E(sigma_new) - E(sigma)
            sigma_old = lattice[i][j];
            sigma_new = sigma_old + delta_dist(gen); 
            
            //E(sigma_new)
            delta_e += (-cos(sigma_new - n0));
            delta_e += (-cos(sigma_new - n1));
            delta_e += (-cos(sigma_new - n2));
            delta_e += (-cos(sigma_new - n3));
            //-E(sigma_old)
            delta_e -= (-cos(sigma_old - n0));
            delta_e -= (-cos(sigma_old - n1));
            delta_e -= (-cos(sigma_old - n2));
            delta_e -= (-cos(sigma_old - n3));


        }
    }

};
inline void over_relaxation_sweep();
inline void cluster_sweep();


///////////////////////
//Standard metropolis//
///////////////////////
void metropolis(int n_steps_therm, int n_steps_prod, int side_length, double beta, std::string outfile, int outfreq, int conf_outfreq){
    //store lattice as 2d array
    std::vector<std::vector<double> > lattice;
    
    //Initialize 128 bits of random state
    std::random_device rd;
    std::array<uint32_t,4> seed;
    seed[0] = rd();
    seed[1] = rd();
    seed[2] = rd();
    seed[3] = rd();
    xoroshiro128plus gen(seed);


    std::uniform_real_distribution<double> angle_dist(-M_PI, M_PI);
    std::uniform_real_distribution<double> delta_dist(-delta/2.0, delta/2.0);
    std::vector<double> energies;
    
    //setup lattice (hot start)
    for (int i=0; i<side_length; i++){
        std::vector<int> line;
        for (int j=0; j<side_length; j++){
            line.push_back(angle_dist(gen));     
        }
        line.shrink_to_fit();
        lattice.push_back(line);
    }
    lattice.shrink_to_fit();

    //thermalize
    for (int n=0; n<n_steps_therm; n++){
        metropolis_sweep();
    }
    //production run
    for (int n=0; n<n_steps_prod; n++){
        metropolis_sweep();
        if (n%outfreq == 0){
            energies.push_back(compute_energy(lattice));
        }
        if (n%conf_outfreq == 0){
            write_configuration(lattice, outfile + ".conf" + std::to_string(n));
        }
    }
};


