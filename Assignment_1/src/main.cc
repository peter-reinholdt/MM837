#include <iostream>
#include "initialize.h"
#include "integrators.h"

int main(int argc, char** argv){
    if(argc < 2){
        std::cout << "Missing input file.\nUsage: FPU.x inputfile\n";
        return 1;
    }
    configuration con;
    initialize(argv[1], con);
    if (con.Settings.integrator == "leapfrog"){
        leapfrog(con.x, con.p, con.Settings.k, con.Settings.nsteps, con.Settings.dt, con.Settings.ifreqout, con.Settings.outfile);
    } else if (con.Settings.integrator == "velocityVerlet"){
        velocityVerlet(con.x, con.p, con.Settings.k, con.Settings.nsteps, con.Settings.dt, con.Settings.ifreqout, con.Settings.outfile);
    } else {
        std::cout << "Error: No integrator specified in settings\n";
    }
    return 0;
}
