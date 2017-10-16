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
        std::vector<double> k;
        std::string outfile;
    void fromFile(std::string filename);
};


class configuration{
    public:
        settings Settings;
        std::vector<double> x;
        std::vector<double> p;
};


void initialize(std::string filename, configuration& con);
