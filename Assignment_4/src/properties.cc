#include <vector>
#include <fstream>
#include <cmath>


inline void to_interval(double& x){
    //0..2*PI interval
    if (x< -1.0*M_PI){
        x += 2*M_PI;
    } else if (x > 1.0*M_PI){
        x -= 2*M_PI;
    }
    
}


void lattice_to_interval(std::vector<std::vector<double> >& lattice){
    int L = lattice.size();
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            to_interval(lattice[i][j]);
        }
    }
}


double compute_energy(std::vector<std::vector<double> >& lattice){
    int L = lattice.size();
    int i_plus, j_plus;
    double energy = 0.0;
    double sigma;
    std::vector<double> nb_angles = {0.0, 0.0};
    /* Doing     or 
     * [ ][0][ ]   [ ][0][ ]
     * [ ][x][0]   [0][x][0]
     * [ ][ ][ ]   [ ][0][ ]
     * is equivalent (after normalization with either 2*L*L  or 4*L*L
     * since the cos() is an even function and the interaction i-j is the same as j-i
     */
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            i_plus  = (i != L-1) ? i+1: 0;
            j_plus  = (j != L-1) ? j+1: 0;
            
            //neighboring angles
            nb_angles[0] = lattice[i_plus][j];
            nb_angles[1] = lattice[i][j_plus];
            //delta_E = E(sigma_new) - E(sigma)
            sigma = lattice[i][j];
            energy -= cos(sigma - nb_angles[0]) + cos(sigma - nb_angles[1]);
        }
    }
    return energy / (2*L*L);
};



double rho(std::vector<double> x, int t){
    double x_mean_0 = 0.0;
    double x_mean_t = 0.0;
    double r = 0.0;
    int n = x.size() - t;
    //get averages
    for (int t0=0; t0<n; t0++){
        x_mean_0 += x[t0];
        x_mean_t += x[t0+t];
    }
    x_mean_0 /= n;
    x_mean_t /= n;
    
    //get t-correlations
    for (int t0=0; t0<n; t0++){
        r += (x[t0] - x_mean_0) * (x[t0+t] - x_mean_t);
    }
    r /= n;
    return r;
};


std::vector<double> compute_autocorr(std::vector<double> x){
    //get vector r(t);
    int t_max = 1000;
    /*
    if ((int)x.size() / 10 < t_max){
        t_max = x.size() / 10;
    }
    */
    std::vector<double> r;

    for (int t=0; t<t_max; t++){
        r.push_back(rho(x,t));
    }
    double norm = 1.0/r[0];
    for (int t=0; t<t_max; t++){
        r[t] *= norm;
    }
    return r;
}


void write_properties(std::vector<double> energies, std::string outfile){
    std::ofstream of;
    auto r = compute_autocorr(energies);
    double tau_int = 0.0;
    of.open(outfile);
    for (int t = 0; t<(int)r.size(); t++){
       tau_int += r[t];
       of << t << " " << r[t] << " " << tau_int << "\n";  
    }
    of.close();
}



void write_energies(std::vector<double> energies, std::string outfile){
    std::ofstream of;
    of.open(outfile);
    for (auto & element: energies){
       of << element << "\n";
    }
    of.close();
}


void write_configuration(std::vector<std::vector<double> >& lattice, std::string outfile){
    std::ofstream of;
    of.open(outfile);
    int L = lattice.size();
    lattice_to_interval(lattice);
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            of << lattice[i][j] << " ";
        }
        of << "\n";
    }
    of.flush();
    of.close();
}
