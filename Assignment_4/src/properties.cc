#include <vector>
#include <fstream>

inline int delta_function(int i, int j){
    if (i==j){
        return 1;
    }else{
        return 0;
    }
}


double compute_energy(std::vector<std::vector<int> >& lattice){
    int L = lattice.size();
    double energy = 2*L*L;
    int i_plus, j_plus;
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            i_plus  = (i != L-1) ? i+1: 0;
            j_plus  = (j != L-1) ? j+1: 0;
            energy -= delta_function(lattice[i][j], lattice[i_plus][j]);
            energy -= delta_function(lattice[i][j], lattice[i][j_plus]);
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
    if ((int)x.size() / 10 < t_max){
        t_max = x.size() / 10;
    }
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
       of << t << ", " << r[t] << ", " << tau_int << "\n";  
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


void write_configuration(std::vector<std::vector<int> >& lattice, std::string outfile){
    std::ofstream of;
    of.open(outfile);
    int L = lattice.size();
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            of << lattice[i][j] << " ";
        }
        of << "\n";
    }
    of.flush();
    of.close();
}
