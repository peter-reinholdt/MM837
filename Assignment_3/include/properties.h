#include <vector>
#include <fstream>



int compute_energy(std::vector<std::vector<int> >& lattice);
void write_properties(std::vector<int> energies, std::string outfile);
