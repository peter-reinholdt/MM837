#include <vector>
#include "toml.h"
#include <random>
#include <iostream>
#include <cmath>
#include "xoroshiro128plus.h"
#include "statistics.h"
#include "properties.h"





inline void metropolis_sweep(std::vector<std::vector<double> >& lattice, 
                             xoroshiro128plus& gen, 
                             std::uniform_real_distribution<double>& delta_dist,
                             std::uniform_real_distribution<double>& real_dist,
                             double beta){
    int L = lattice.size();
    int i_plus, j_plus, i_minus, j_minus;
    std::vector<double> nb_angles = {0.0, 0.0, 0.0, 0.0};
    double sigma_old, sigma_new;
    double delta_E, p_accept;
    N_ATTEMPTED_FLIPS += L*L;

    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            delta_E = 0.0;
            i_plus  = (i != L-1) ? i+1: 0;
            i_minus = (i != 0  ) ? i-1: L-1;
            j_plus  = (j != L-1) ? j+1: 0;
            j_minus = (j != 0  ) ? j-1: L-1;
            
            //neighboring angles
            nb_angles[0] = lattice[i_plus][j];
            nb_angles[1] = lattice[i_minus][j];
            nb_angles[2] = lattice[i][j_plus];
            nb_angles[3] = lattice[i][j_minus];

            //delta_E = E(sigma_new) - E(sigma)
            sigma_old = lattice[i][j];
            sigma_new = sigma_old - delta_dist(gen); 
            
            for (int nb=0; nb<4; nb++){
                delta_E += cos(sigma_old - nb_angles[nb]) - cos(sigma_new - nb_angles[nb]);
            }
            
            if (delta_E <= 0){
                //good change
                lattice[i][j] = sigma_new;
                N_ACCEPTED_FLIPS++;
            }else{
                p_accept = exp(-beta*delta_E);
                if (p_accept > real_dist(gen)){
                    lattice[i][j] = sigma_new;
                    N_ACCEPTED_FLIPS++;
                }
            }
        }
    }
};



inline void microcanonical_sweep(std::vector<std::vector<double> >& lattice){
    int L = lattice.size();
    int i_plus, j_plus, i_minus, j_minus;
    int n0, n1, n2, n3;

    std::vector<double> Vx;

    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){

            i_plus  = (i != L-1) ? i+1: 0;
            i_minus = (i != 0  ) ? i-1: L-1;
            j_plus  = (j != L-1) ? j+1: 0;
            j_minus = (j != 0  ) ? j-1: L-1;
            
            //neighboring angles
            
        }
    }


};
inline void cluster_sweep(){};


///////////////////////
//Standard metropolis//
///////////////////////
void metropolis(int n_steps_therm, int n_steps_prod, int side_length, double beta, std::string outfile, int outfreq, int conf_outfreq, double delta){
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


    std::uniform_real_distribution<double> angle_dist(0.0, 2*M_PI);
    std::uniform_real_distribution<double> delta_dist(-delta/2.0, delta/2.0);
    std::uniform_real_distribution<double> real_dist(0.0, 1.0);
    std::vector<double> energies;
    
    //setup lattice (hot start)
    for (int i=0; i<side_length; i++){
        std::vector<double> line;
        for (int j=0; j<side_length; j++){
            line.push_back(angle_dist(gen));     
        }
        line.shrink_to_fit();
        lattice.push_back(line);
    }
    lattice.shrink_to_fit();
    

    
    //thermalize
    std::cout << lattice[0][0] << std::endl;
    for (int n=0; n<n_steps_therm; n++){
        metropolis_sweep(lattice, gen, delta_dist, real_dist, beta);

        std::cout << lattice[0][0] << std::endl;

    }
    //production run
    for (int n=0; n<n_steps_prod; n++){
        metropolis_sweep(lattice, gen, delta_dist, real_dist, beta);
        if (n%outfreq == 0){
        }
        if (n%conf_outfreq == 0){
            write_configuration(lattice, outfile + ".conf" + std::to_string(n));
        }
    }
};


