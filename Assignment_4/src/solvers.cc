#include <vector>
#include "toml.h"
#include <random>
#include <iostream>
#include <cmath>
#include "xoroshiro128plus.h"
#include "statistics.h"
#include "properties.h"
#include <stack>
#include <set>




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
    std::vector<double> nb_angles = {0.0, 0.0, 0.0, 0.0};
    double Vxx, Vxy, Vx_sq;
    double x, y, x_new, y_new;
    double sigma;

    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            Vxx = 0.0;
            Vxy = 0.0;
            sigma = lattice[i][j];
            x = cos(sigma);
            y = sin(sigma);

            i_plus  = (i != L-1) ? i+1: 0;
            i_minus = (i != 0  ) ? i-1: L-1;
            j_plus  = (j != L-1) ? j+1: 0;
            j_minus = (j != 0  ) ? j-1: L-1;
            
            //neighboring angles
            nb_angles[0] = lattice[i_plus][j];
            nb_angles[1] = lattice[i_minus][j];
            nb_angles[2] = lattice[i][j_plus];
            nb_angles[3] = lattice[i][j_minus];

            for (int nb=0; nb<4; nb++){
                Vxx += cos(nb_angles[nb]);
                Vxy += sin(nb_angles[nb]);
            }
            Vx_sq = Vxx*Vxx + Vxy*Vxy; 
            x_new = 2 * (x * Vxx + y * Vxy) / Vx_sq * Vxx - x;
            y_new = 2 * (x * Vxx + y * Vxy) / Vx_sq * Vxy - y;
            lattice[i][j] = atan2(y_new, x_new);
        }
    }

};


inline void cluster_sweep(std::vector<std::vector<double> >& lattice,
                          xoroshiro128plus& gen,
                          std::uniform_real_distribution<double>& angle_dist,
                          std::uniform_real_distribution<double>& real_dist,
                          std::uniform_int_distribution<int>& L_dist,
                          double beta){
    int i_plus, j_plus, i_minus, j_minus;
    int L = lattice.size();
    double p_add, sx, sy;
    N_ATTEMPTED_FLIPS += L*L;

    //make stack of neighbor spins
    std::stack<std::pair<int,int> > cluster_buffer;
    std::pair<int,int> current_site;

    //make set of sites currently in the cluster
    std::set<std::pair<int,int> > cluster_set;
    
    //pick random vector r
    double theta = angle_dist(gen);
    double rx = cos(theta);
    double ry = sin(theta);

    //pick starting seed
    int i = L_dist(gen);
    int j = L_dist(gen);

    //add seed site to accepted spins
    N_ACCEPTED_FLIPS++;
    cluster_set.insert(std::make_pair(i,j));
    cluster_buffer.push(std::make_pair(i, j));
    

    while(!cluster_buffer.empty()){
        std::vector<std::pair<int,int> > neighbors;
        //look at the top of the stack, and find the neighbors
        current_site = cluster_buffer.top();
        cluster_buffer.pop();

        i = current_site.first; 
        j = current_site.second; 
        sx = cos(lattice[i][j])*rx + sin(lattice[i][j])*ry;

        i_plus  = (i != L-1) ? i+1: 0;
        i_minus = (i != 0  ) ? i-1: L-1;
        j_plus  = (j != L-1) ? j+1: 0;
        j_minus = (j != 0  ) ? j-1: L-1;

        neighbors.push_back(std::make_pair(i_minus, j));
        neighbors.push_back(std::make_pair(i_plus, j));
        neighbors.push_back(std::make_pair(i, j_minus));
        neighbors.push_back(std::make_pair(i, j_plus));

        for (int n=0; n<(int)neighbors.size(); n++){
            i = neighbors[n].first;
            j = neighbors[n].second;
            if (cluster_set.count(neighbors[n]) == 0){
                // neighbor is *not* in cluster already
                // check if it should be added
                sy = cos(lattice[i][j]*rx) + sin(lattice[i][j])*ry;
                p_add = 1.0 - exp(-2*beta*sx*sy); 
                if (real_dist(gen) < p_add){
                    N_ACCEPTED_FLIPS++;
                    cluster_set.insert(neighbors[n]);
                    cluster_buffer.push(neighbors[n]);
                }
           }
        }
    }
    for (auto site : cluster_set){
        //flip all spins in the cluster
        i = site.first;
        j = site.second;
        sx = cos(lattice[i][j])*rx + sin(lattice[i][j])*ry;
        lattice[i][j] = atan2(lattice[i][j] - 2*sx*ry, lattice[i][j] - 2*sx*rx);
    }
}
        




////////////////////////////////////////
//Standard metropolis + microcanonical//
////////////////////////////////////////
void metropolis(int n_steps_therm, int n_steps_prod, int side_length, double beta, std::string outfile, int outfreq, int conf_outfreq, double delta, int n_mc_sweep){
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


    std::uniform_real_distribution<double> angle_dist(0.0, 2.0*M_PI);
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
    for (int n=0; n<n_steps_therm; n++){
        metropolis_sweep(lattice, gen, delta_dist, real_dist, beta);
        //do n_mc_sweep of microcanonical sweeps
        for (int mc=0; mc<n_mc_sweep; mc++){microcanonical_sweep(lattice);}
        if (n%100 == 0){
            std::cout << "\rThermalizing, " << std::setprecision(3) << 100 * (double)n / (double)n_steps_therm << "%" << std::flush;
        }
    }
    std::cout << "\rThermalizing, ...Done!";
    std::cout << std::endl;
    //production run
    for (int n=0; n<n_steps_prod; n++){
        metropolis_sweep(lattice, gen, delta_dist, real_dist, beta);
        for (int mc=0; mc<n_mc_sweep; mc++){microcanonical_sweep(lattice);}
        if (n%outfreq == 0){
            energies.push_back(compute_energy(lattice));
        }
        if (n%conf_outfreq == 0){
            write_configuration(lattice, outfile + ".conf" + std::to_string(n));
        }
        if (n%100 == 0){
            std::cout << "\rRunning production, " << std::setprecision(3) << 100 * (double)n / (double)n_steps_prod << "%" << std::flush;
        }
    }
    std::cout << "\rRunning production, ...Done!";
    std::cout << std::endl;

    std::cout << "Beta: " << beta << " Acceptance ratio: " << (double) N_ACCEPTED_FLIPS / (double) N_ATTEMPTED_FLIPS << std::endl;
    //write properties and final configuration to file
    write_properties(energies, outfile + ".autocorr");
    write_energies(energies, outfile + ".energies");
    write_configuration(lattice, outfile + ".conf");

};

/////////////////////
//Cluster algorithm//
/////////////////////
void cluster(int n_steps_therm, int n_steps_prod, int side_length, double beta, std::string outfile, int outfreq, int conf_outfreq, int n_mc_sweep){
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

    std::uniform_real_distribution<double> angle_dist(0.0, 2.0*M_PI);
    std::uniform_real_distribution<double> real_dist(0.0, 1.0);
    std::uniform_int_distribution<int> L_dist(0, side_length-1);
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
    for (int n=0; n<n_steps_therm; n++){



        cluster_sweep(lattice, gen, angle_dist, real_dist, L_dist, beta);
        //do n_mc_sweep of microcanonical sweeps
        for (int mc=0; mc<n_mc_sweep; mc++){microcanonical_sweep(lattice);}
        if (n%100 == 0){
            std::cout << "\rThermalizing, " << std::setprecision(3) << 100 * (double)n / (double)n_steps_therm << "%" << std::flush;
        }
    }
    std::cout << "\rThermalizing, ...Done!";
    std::cout << std::endl;
    //production run
    for (int n=0; n<n_steps_prod; n++){
        cluster_sweep(lattice, gen, angle_dist, real_dist, L_dist, beta);
        for (int mc=0; mc<n_mc_sweep; mc++){microcanonical_sweep(lattice);}
        if (n%outfreq == 0){
            energies.push_back(compute_energy(lattice));
        }
        if (n%conf_outfreq == 0){
            write_configuration(lattice, outfile + ".conf" + std::to_string(n));
        }
        if (n%100 == 0){
            std::cout << "\rRunning production, " << std::setprecision(3) << 100 * (double)n / (double)n_steps_prod << "%" << std::flush;
        }
    }
    std::cout << "\rRunning production, ...Done!";
    std::cout << std::endl;

    std::cout << "Beta: " << beta << " Acceptance ratio: " << (double) N_ACCEPTED_FLIPS / (double) N_ATTEMPTED_FLIPS << std::endl;
    //write properties and final configuration to file
    write_properties(energies, outfile + ".autocorr");
    write_energies(energies, outfile + ".energies");
    write_configuration(lattice, outfile + ".conf");

};
