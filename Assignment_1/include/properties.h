#include <vector>
#include <fstream>
#include "initialize.h"


double computeMeanPsq(const std::vector<double>& p);
double computeKineticEnergy(const std::vector<double>& p);
double computePotentialEnergy(const std::vector<double>& x, const std::vector<double> k);
void writeProperties(const std::vector<double>& x, const std::vector<double>& p, const std::vector<double>& k, FILE * outfile);
void dumpCoordinates(const std::vector<double>& x, const std::vector<double>& p, FILE * coords);
