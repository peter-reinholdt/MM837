#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include "forces.h"


void leapfrog(std::vector<double>& x, std::vector<double>& p, int nsteps, double dt, int ifreqout, std::string outfile){
    //kick-drift-kick form
    int N = x.size();
    std::vector<double> forces;
    forces.resize(N);
    computeForces(x, forces);
    for(int n=0; n<nsteps; n++){
    /*    if( n % ifreqout == 0){
           // format: timestep,x0,x1,...p0,p1
            std::cout << n << "\n";
        }*/
        for(int i=0; i<N; i++){
            //p_{i+0.5} = p_i + a_i * 0.5 * dt
            p[i] += forces[i] * 0.5 * dt;
            //x_{i+1}   = x_i + p_{i+0.5} * dt
            x[i] += p[i] * dt;
        }
        computeForces(x, forces);
        for(int i=0; i<N; i++){
            //p_{i+1}   = p_{i+0.5} + a_{i+1} * 0.5 * dt
            p[i] += forces[i] * 0.5 * dt;
        }
    }
}
