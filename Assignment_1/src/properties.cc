#include <vector>
#include <fstream>
#include <iostream>
#include "initialize.h"


double computeMeanPsq(const std::vector<double>& p){
    double meanPsq = 0.0;
    for(int i=0; i<p.size(); i++){
        meanPsq += p[i]*p[i];
    } 
    return meanPsq / p.size();
}


double computeKineticEnergy(const std::vector<double>& p){
    double Ekin = 0.0;
    for(int i=0; i<p.size(); i++){
        Ekin += 0.5 * p[i] * p[i];
    }
    return Ekin;
}


double computePotentialEnergy(const std::vector<double>& x, const std::vector<double> k){
    double Epot = 0.0;
    double x1;
    double x2;
    double x3;
    double x4;
    //potential from i=1..N-1
    for(int i=1; i<x.size(); i++){
        x1 = x[i] - x[i-1];
        x2 = x1 * x1;
        x3 = x2 * x1; 
        x4 = x3 * x1;
        Epot += k[0] / 2.0 * x2 + k[1] / 3.0 * x3 + k[2] / 4.0 * x4;
    }
    //potential from i=0
    x1 = x[0] - x[x.size()-1];
    x2 = x1 * x1;
    x3 = x2 * x1; 
    x4 = x3 * x1;
    Epot += k[0] / 2.0 * x2 + k[1] / 3.0 * x3 + k[2] / 4.0 * x4;
    return Epot;
}


void writeProperties(const std::vector<double>& x, const std::vector<double>& p, const std::vector<double>& k, FILE * outfile){
    double Ekin = computeKineticEnergy(p);
    double Epot = computePotentialEnergy(x, k);
    double meanPsq = computeMeanPsq(p);
                  //<p^2>            E_tot                 Ekin           Epot
    fprintf(outfile, "%12.6f %12.6f %12.6f %12.6f\n", meanPsq, Ekin+Epot, Ekin, Epot);
}
