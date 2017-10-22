#include <vector>
#include <fstream>
#include <string>
#include "forces.h"
#include "properties.h"

void leapfrog(std::vector<double>& x, std::vector<double>& p, const std::vector<double> k, const int nsteps, const double dt, const int ifreqout, const std::string outfile){
    //kick-drift-kick form
    int N = x.size();
    std::vector<double> forces;
    forces.resize(N);
    computeForces(x, forces, k);
    
    FILE * of;    
    of = fopen(outfile.c_str(), "w");

    FILE * cf;    
    cf = fopen("coords.dat", "w");

    for(int n=0; n<nsteps; n++){
        if (n%ifreqout == 0){
            writeProperties(x, p, k, of);
            dumpCoordinates(x, p, cf);
        }
        for(int i=0; i<N; i++){
            //p_{i+0.5} = p_i + a_i * 0.5 * dt
            p[i] += forces[i] * 0.5 * dt;
            //x_{i+1}   = x_i + p_{i+0.5} * dt
            x[i] += p[i] * dt;
        }
        computeForces(x, forces, k);
        for(int i=0; i<N; i++){
            //p_{i+1}   = p_{i+0.5} + a_{i+1} * 0.5 * dt
            p[i] += forces[i] * 0.5 * dt;
        }
    }
    fclose(of);
    fclose(cf);
}



void velocityVerlet(std::vector<double>& x, std::vector<double>& p, const std::vector<double> k, const int nsteps, const double dt, const int ifreqout, const std::string outfile){
    int N = x.size();
    std::vector<double> forces;
    forces.resize(N);
    std::vector<double> oldforces;
    oldforces.resize(N);
    computeForces(x, forces, k);

    FILE * of;    
    of = fopen(outfile.c_str(), "w");

    FILE * cf;    
    cf = fopen("coords.dat", "w");
    
    
    for(int n=0; n<nsteps; n++){
        if (n%ifreqout == 0){
            writeProperties(x, p, k, of);
            dumpCoordinates(x, p, cf);
        }
        for(int i=0; i<N; i++){
            x[i] += p[i] * dt + 0.5 * forces[i] * dt * dt;
        }
        oldforces = forces;
        computeForces(x, forces, k);
        for(int i=0; i<N; i++){
            p[i] += 0.5 * (forces[i] + oldforces[i] ) * dt;
        }
    }
    fclose(of);
    fclose(cf);
}


