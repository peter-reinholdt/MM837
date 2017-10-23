#include <iostream>
#include <chrono>
#include "initialize.h"
#include "integrators.h"


int main(int argc, char** argv){
    if(argc < 2){
        std::cout << "Missing input file.\nUsage: spring [inputfile]\n";
        return 1;
    }
    
    //initialize start so compiler is happy
    auto start = std::chrono::system_clock::now();
    
    configuration con;
    initialize(argv[1], con);
    if (con.Settings.integrator == "leapfrog"){
        //real start now
        start = std::chrono::system_clock::now();
        leapfrog(con.x, con.p, con.Settings.k, con.Settings.nsteps, con.Settings.dt, con.Settings.ifreqout, con.Settings.outfile);
    } else if (con.Settings.integrator == "velocityVerlet"){
        //or now
        start = std::chrono::system_clock::now();
        velocityVerlet(con.x, con.p, con.Settings.k, con.Settings.nsteps, con.Settings.dt, con.Settings.ifreqout, con.Settings.outfile);
    } else {
        std::cout << "Error: No integrator specified in settings\n";
        return 1;
    }
    
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "N = "<< con.Settings.Nparticles << " Elapsed time (s): " << elapsed_seconds.count() << "\n";
    return 0;
}
