#include <vector>
#include <stack>
#include <time.h>
#include "toml.h"
#include <random>
#include <iostream>
#include "properties.h"
#include "statistics.h"
#include "xoroshiro128plus.h"
#include <array>


inline int delta_function(int i, int j){
    if (i==j){
        return 1;
    }else{
        return 0;
    }
}

inline void metropolis_sweep(std::vector<std::vector<int> >& lattice,
                             double beta,
                             int q_max,
                             xoroshiro128plus& gen,
                             std::uniform_int_distribution<int>& q_dist, 
                             std::uniform_real_distribution<double>& double_dist){
    int delta_E, sigma_old, sigma_new;
    int L = lattice.size();
    int i_plus, i_minus, j_plus, j_minus;
    std::vector<double> exp_minus_beta_E;
    double p_accept;

    N_ATTEMPTED_FLIPS += L*L;
    //(local) delta_E on spin flip is -4, -3, -2, -1, 0, 1, 2, 3, 4.
    //we always accept delta_E <= 0; for 0, 1, 2, 3, 4
    //make a list that we can index with delta_E to get p_acc
    //this avoids ~(L**2) calculations of exp(-beta*i)
    for (int i=0; i<5; i++){
       exp_minus_beta_E.push_back(exp(-beta*i)); 
    }

    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            
            i_plus  = (i != L-1) ? i+1: 0;
            i_minus = (i != 0  ) ? i-1: L-1;
            j_plus  = (j != L-1) ? j+1: 0;
            j_minus = (j != 0  ) ? j-1: L-1;
            

            sigma_old = lattice[i][j];
            sigma_new = q_dist(gen);
            //
            //NOTE: q_dist gives 0,1,...,q-2
            // in total (q-1) values
            // if sigma_new is drawn to be the current value  
            // instead set it to q-1
            // so we map to the same values
            // but don't do stupid work in 
            // drawing extra random numbers
            if (sigma_new == sigma_old){
                sigma_new=q_max-1;
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
                N_ACCEPTED_FLIPS++;
            } else {
                p_accept = exp_minus_beta_E[delta_E]; //exp(-beta*delta_E);
                if (p_accept > double_dist(gen)){
                    //accept with p_accept probability otherwise
                    lattice[i][j] = sigma_new;
                    N_ACCEPTED_FLIPS++;
                }
            }
        }
    }
};


inline void cluster_sweep(std::vector<std::vector<int> >& lattice,
                          double beta,
                          xoroshiro128plus& gen,
                          std::uniform_int_distribution<int>& q_dist,
                          std::uniform_int_distribution<int>& L_dist,
                          std::uniform_real_distribution<double>& double_dist){
    int i_plus, j_plus, i_minus, j_minus;
    int L = lattice.size();
    N_ATTEMPTED_FLIPS += L*L;
    int sigma_old, sigma_new;
    std::stack<std::pair<int,int> > cluster_buffer;
    std::pair<int,int> current_site;
    double p_add = 1 - exp(-beta);
    
    //pick starting seed
    int i = L_dist(gen);
    int j = L_dist(gen);

    //pick *new* value of q
    sigma_old = lattice[i][j];
    do{
        sigma_new = q_dist(gen);
    } while (sigma_new == sigma_old);

    //flip spin of original and add to stack
    lattice[i][j] = sigma_new;
    N_ACCEPTED_FLIPS++;
    cluster_buffer.push(std::make_pair(i, j));

    while(!cluster_buffer.empty()){
        std::vector<std::pair<int,int> > neighbours;
        //look at the top of the stack, and find the neighbours
        current_site = cluster_buffer.top();
        cluster_buffer.pop();

        i = current_site.first; 
        j = current_site.second; 

        i_plus  = (i != L-1) ? i+1: 0;
        i_minus = (i != 0  ) ? i-1: L-1;
        j_plus  = (j != L-1) ? j+1: 0;
        j_minus = (j != 0  ) ? j-1: L-1;

        neighbours.push_back(std::make_pair(i_minus, j));
        neighbours.push_back(std::make_pair(i_plus, j));
        neighbours.push_back(std::make_pair(i, j_minus));
        neighbours.push_back(std::make_pair(i, j_plus));
        for (int n=0; n<(int)neighbours.size(); n++){
            i = neighbours[n].first;
            j = neighbours[n].second;
            if (lattice[i][j] == sigma_old){
                if (double_dist(gen) < p_add){
                    N_ACCEPTED_FLIPS++;
                    //flip spin now, so it is not considered later
                    lattice[i][j] = sigma_new;
                    cluster_buffer.push(std::make_pair(i,j));
                }
            }
        }
    }

};



///////////////////////
//Standard metropolis//
///////////////////////
void metropolis(int n_steps_therm, int n_steps_prod, int side_length, int potts_q, double beta, std::string outfile, int outfreq, int conf_outfreq){
    //store lattice as 2d array
    std::vector<std::vector<int> > lattice;
    //std::mt19937 gen(time(NULL));
    //
    std::random_device rd;
    std::array<uint32_t,4> seed;
    seed[0] = rd();
    seed[1] = rd();
    seed[2] = rd();
    seed[3] = rd();
    xoroshiro128plus gen(seed);


    std::uniform_int_distribution<int> q_dist(0, potts_q-2);
    std::uniform_int_distribution<int> q_dist_start(0, potts_q-1);
    std::uniform_real_distribution<double> real_dist(0.0, 1.0);
    std::vector<double> energies;
    
    //setup lattice (hot start)
    for (int i=0; i<side_length; i++){
        std::vector<int> line;
        for (int j=0; j<side_length; j++){
            line.push_back(q_dist_start(gen));     
        }
	line.shrink_to_fit();
        lattice.push_back(line);
    }
    lattice.shrink_to_fit();


    //thermalize
    for (int n=0; n<n_steps_therm; n++){
        metropolis_sweep(lattice, beta, potts_q, gen, q_dist, real_dist);
    }
    //production run
    for (int n=0; n<n_steps_prod; n++){
        metropolis_sweep(lattice, beta, potts_q, gen, q_dist, real_dist);
        if (n%outfreq == 0){
            energies.push_back(compute_energy(lattice));
        }
        if (n%conf_outfreq == 0){
            write_configuration(lattice, outfile + ".conf" + std::to_string(n));
        }
    }

    std::cout << "Beta: " << beta << " Acceptance ratio: " << (double) N_ACCEPTED_FLIPS / (double) N_ATTEMPTED_FLIPS << std::endl;

    //write properties and final configuration to file
    write_properties(energies, outfile + ".autocorr");
    write_energies(energies, outfile + ".energies");
    write_configuration(lattice, outfile + ".conf");
};


////////////////////
//Standard cluster//
////////////////////
void cluster(int n_steps_therm, int n_steps_prod, int side_length, int potts_q, double beta, std::string outfile, int outfreq, int conf_outfreq){
    //store lattice as 2d array
    std::vector<std::vector<int> > lattice;
    std::random_device rd;
    std::array<uint32_t,4> seed;
    seed[0] = rd();
    seed[1] = rd();
    seed[2] = rd();
    seed[3] = rd();
    xoroshiro128plus gen(seed);
    std::uniform_int_distribution<int> q_dist(0, potts_q-1);
    std::uniform_int_distribution<int> L_dist(0, side_length-1);
    std::uniform_real_distribution<double> real_dist(0.0, 1.0);
    std::vector<double> energies;
    
    //setup lattice (hot start)
    for (int i=0; i<side_length; i++){
        std::vector<int> line;
        for (int j=0; j<side_length; j++){
            line.push_back(q_dist(gen));     
        }
        lattice.push_back(line);
    }


    //thermalize
    
    for (int n=0; n<n_steps_therm; n++){
        cluster_sweep(lattice, beta, gen, q_dist, L_dist, real_dist);
    }
    //production run
    for (int n=0; n<n_steps_prod; n++){
        cluster_sweep(lattice, beta, gen, q_dist, L_dist, real_dist);
        if (n%outfreq == 0){
            energies.push_back(compute_energy(lattice));
        }
        if (n%conf_outfreq == 0){
            write_configuration(lattice, outfile + ".conf" + std::to_string(n));
        }
    }

    std::cout << "Beta: " << beta << " Acceptance ratio: " << (double) N_ACCEPTED_FLIPS / (double) N_ATTEMPTED_FLIPS << std::endl;


    //write properties and final configuration to file
    write_properties(energies, outfile + ".autocorr");
    write_energies(energies, outfile + ".energies");
    write_configuration(lattice, outfile + ".conf");
};



/////////////////////////////
//Hybrid cluster-metropolis//
/////////////////////////////
void hybrid(int n_steps_therm, int n_steps_prod, int side_length, int potts_q, double beta, std::string outfile, int outfreq, int conf_outfreq, int metropolis_freq){
    //store lattice as 2d array
    std::vector<std::vector<int> > lattice;
    //std::mt19937 gen(time(NULL));
    //
    std::random_device rd;
    std::array<uint32_t,4> seed;
    seed[0] = rd();
    seed[1] = rd();
    seed[2] = rd();
    seed[3] = rd();
    xoroshiro128plus gen(seed);


    std::uniform_int_distribution<int> q_dist(0, potts_q-2);
    std::uniform_int_distribution<int> q_dist_start(0, potts_q-1);
    std::uniform_real_distribution<double> real_dist(0.0, 1.0);
    std::uniform_int_distribution<int> L_dist(0, side_length-1);
    std::vector<double> energies;
    
    //setup lattice (hot start)
    for (int i=0; i<side_length; i++){
        std::vector<int> line;
        for (int j=0; j<side_length; j++){
            line.push_back(q_dist_start(gen));     
        }
	line.shrink_to_fit();
        lattice.push_back(line);
    }
    lattice.shrink_to_fit();


    //thermalize
    for (int n=0; n<n_steps_therm; n++){
	if (n%(metropolis_freq+1)==0){
            metropolis_sweep(lattice, beta, potts_q, gen, q_dist, real_dist);
	}else{
            cluster_sweep(lattice, beta, gen, q_dist_start, L_dist, real_dist);
	}
    }
    //production run
    for (int n=0; n<n_steps_prod; n++){
	if (n%(metropolis_freq+1)==0){
            metropolis_sweep(lattice, beta, potts_q, gen, q_dist, real_dist);
	}else{
            cluster_sweep(lattice, beta, gen, q_dist_start, L_dist, real_dist);
	}
        if (n%outfreq == 0){
            energies.push_back(compute_energy(lattice));
        }
        if (n%conf_outfreq == 0){
            write_configuration(lattice, outfile + ".conf" + std::to_string(n));
        }
    }

    std::cout << "Beta: " << beta << " Acceptance ratio: " << (double) N_ACCEPTED_FLIPS / (double) N_ATTEMPTED_FLIPS << std::endl;

    //write properties and final configuration to file
    write_properties(energies, outfile + ".autocorr");
    write_energies(energies, outfile + ".energies");
    write_configuration(lattice, outfile + ".conf");
};
