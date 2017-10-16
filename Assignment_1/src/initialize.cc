#include <string>
#include <iostream>
#include <fstream>
#include "toml.h"


class settings{
    private:
    public:
        int Nparticles;
        int nsteps;
        int ifreqout;
        double delta;
        double dt;
        double m;
        std::vector<double> k;
        std::string outfile;
    void fromFile(std::string filename){
        std::ifstream ifs(filename);
        toml::ParseResult pr = toml::parse(ifs);
        toml::Value v = pr.value;
        //get integrator settings
        const toml::Value& integrators = v.get<toml::Array>("integrators")[0];
        dt = integrators.get<double>("dt");
        nsteps = integrators.get<int>("nsteps");
        //get system details
        const toml::Value& system = v.get<toml::Array>("system")[0];
        Nparticles = system.get<int>("Nparticles");
        delta = system.get<double>("delta");
        k.resize(3);
        k[0] = system.get<double>("k1");
        k[1] = system.get<double>("k2");
        k[2] = system.get<double>("k3");
        
        //get output details
        const toml::Value& output = v.get<toml::Array>("output")[0];
        ifreqout = output.get<int>("ifreqout");
        outfile = output.get<std::string>("outfile");
        ifs.close();
    }
};

class configuration{
    public:
        settings Settings;
        std::vector<double> x;
        std::vector<double> p;
};


void initialize(std::string filename, configuration& con){
    //read input file
    settings st;
    st.fromFile(filename);
    con.Settings = st;
    con.x.resize(st.Nparticles);
    con.p.resize(st.Nparticles);
    //set initial p, x
    for(int i=0; i<st.Nparticles; i++){
        con.x[i] = 0.0;
        con.p[i] = sin((2*M_PI*(i+st.delta)/st.Nparticles));
    }
}
