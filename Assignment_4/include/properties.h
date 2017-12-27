#include <vector>
#include <fstream>



double compute_energy(std::vector<std::vector<double> >& lattice);
void write_properties(std::vector<double> energies, std::string outfile);
void write_energies(std::vector<double> energies, std::string outfile);
void write_configuration(std::vector<std::vector<double> >& lattice, std::string outfile);
