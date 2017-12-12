#include <vector>
#include <fstream>



double compute_energy(std::vector<std::vector<int> >& lattice);
void write_properties(std::vector<double> energies, std::string outfile);
void write_energies(std::vector<double> energies, std::string outfile);
void write_configuration(std::vector<std::vector<int> >& lattice, std::string outfile);
