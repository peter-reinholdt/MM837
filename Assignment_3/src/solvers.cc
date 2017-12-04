#include <vector>
#include <time.h>
#include "toml.h"
#include <random>
#include <iostream>
#include "properties.h"

inline int delta_function(int i, int j){
    if (i==j){
        return 1;
    }else{
        return 0;
    }
}

inline void metropolis_sweep(std::vector<std::vector<int> >& lattice,
                             double beta,
                             std::mt19937& gen,
                             std::uniform_int_distribution<int>& int_dist, 
                             std::uniform_real_distribution<double>& double_dist){
    int delta_E, sigma_old, sigma_new;
    int L = lattice.size();
    int i_plus, i_minus, j_plus, j_minus;
    double p_accept;

    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            i_plus  = (i != L-1) ? i+1: 0;
            i_minus = (i != 0  ) ? i-1: L-1;
            j_plus  = (j != L-1) ? j+1: 0;
            j_minus = (j != 0  ) ? j-1: L-1;
            
            sigma_old = lattice[i][j];
            sigma_new = int_dist(gen);
            while (sigma_old == sigma_new){
                sigma_new = int_dist(gen);
            }
            //change in energy is delta_E = E(sigma') - E(sigma)
            //local energy is sum_{i,j} (1 - delta(i,j))
            //but we want differences, 
            //so in the change we don't to include the 1.
            delta_E = 0;
            delta_E += ( - delta_function(sigma_new, lattice[i_plus][j]));
            delta_E += ( - delta_function(sigma_new, lattice[i_minus][j]));
            delta_E += ( - delta_function(sigma_new, lattice[i][j_plus]));
            delta_E += ( - delta_function(sigma_new, lattice[i][j_minus]));
            // plus energy after flip
            delta_E -= ( - delta_function(sigma_old, lattice[i_plus][j]));
            delta_E -= ( - delta_function(sigma_old, lattice[i_minus][j]));
            delta_E -= ( - delta_function(sigma_old, lattice[i][j_plus]));
            delta_E -= ( - delta_function(sigma_old, lattice[i][j_minus]));
            
            //check for acceptance
            if (delta_E <= 0){
                //accept if energy is lower
                lattice[i][j] = sigma_new;
            } else {
                p_accept = exp(-beta*delta_E);
                if (p_accept > double_dist(gen)){
                    //accept with p_accept probability otherwise
                    lattice[i][j] = sigma_new;
                }
            }
        }
    }
};


void metropolis(int n_steps_therm, int n_steps_prod, int side_length, int potts_q, double beta, std::string outfile, int outfreq){
    //store lattice as 2d array
    std::vector<std::vector<int> > lattice;
    std::mt19937 gen(time(NULL));
    std::uniform_int_distribution<int> int_dist(0, potts_q-1);
    std::uniform_real_distribution<double> real_dist(0.0, 1.0);
    std::vector<int> energies;
    
    //setup lattice (hot start)
    for (int i=0; i<side_length; i++){
        std::vector<int> line;
        for (int j=0; j<side_length; j++){
            line.push_back(int_dist(gen));     
        }
        lattice.push_back(line);
    }


    //thermalize
    for (int n=0; n<n_steps_therm; n++){
        metropolis_sweep(lattice, beta, gen, int_dist, real_dist);
    }
    //production run
    for (int n=0; n<n_steps_prod; n++){
        metropolis_sweep(lattice, beta, gen, int_dist, real_dist);
        if (n%outfreq == 0){
            energies.push_back(compute_energy(lattice));
        }
    }
    //write properties to file
    write_energy(energies, outfile);
    /* 
    for (int i=0; i<side_length; i++){
        for (int j=0; j<side_length; j++){
            std::cout << lattice[i][j];
        }
        std::cout << "\n";
    }
    */
};



void cluster(int n_steps_therm, int n_steps_prod, int side_length, int potts_q, double beta, std::string outfile, int outfreq){};
