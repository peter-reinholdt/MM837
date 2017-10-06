#include <string>
#include <iostream>
#include <fstream>
#include "toml.h"

class settings{
    public:
        int Nparticles;
        int nsteps;
        int ifreqout;
        double delta;
        double dt;
        double m;
        double k1;
        double k2;
        double k3;
        std::string outfile;
    void fromFile(std::string filename);
};


struct configuration{
    settings Settings;
    std::vector<double> x;
    std::vector<double> p;
};


void initialize(std::string filename, configuration& con);
