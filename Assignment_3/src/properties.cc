#include <vector>
#include <fstream>

inline int delta_function(int i, int j){
    if (i==j){
        return 1;
    }else{
        return 0;
    }
}


int compute_energy(std::vector<std::vector<int> >& lattice){
    int L = lattice.size();
    int energy = -(L*L);
    int i_plus, j_plus;
    for (int i=0; i<L-1; i++){
        for (int j=0; j<L+1; j++){
            i_plus  = (i != L-1) ? i+1: 0;
            j_plus  = (j != L-1) ? j+1: 0;
            energy -= delta_function(lattice[i][j], lattice[i_plus][j]);
            energy -= delta_function(lattice[i][j], lattice[i][j_plus]);
        }
    }
    return energy;
};

void write_energy(std::vector<int> energies, std::string outfile){
    std::ofstream of;
    of.open(outfile);
    for (auto & element: energies){
       of << element << "\n"; 
    }
    of.close();
}

