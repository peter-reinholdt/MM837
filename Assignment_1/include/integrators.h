#include<vector>
#include<string>

void leapfrog(std::vector<double>& x, std::vector<double>& p, const std::vector<double> k, int nsteps, double dt, int ifreqout, std::string outfile);
void velocityVerlet(std::vector<double>& x, std::vector<double>& p, const std::vector<double> k, int nsteps, double dt, int ifreqout, std::string outfile);
