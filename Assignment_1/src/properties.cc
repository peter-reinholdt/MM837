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
    //potential from i=0..N-1
    for(int i=0; i<x.size()-1; i++){
        x1 = x[i] - x[i+1];
        x2 = x1 * x1;
        x3 = x2 * x1; 
        x4 = x3 * x1;
        Epot += k[0] / 2.0 * x2 + k[1] / 3.0 * x3 + k[2] / 4.0 * x4;
    }
    //potential from i=0
    x1 = x[x.size()-1] - x[0];
    x2 = x1 * x1;
    x3 = x2 * x1; 
    x4 = x3 * x1;
    Epot += k[0] / 2.0 * x2 + k[1] / 3.0 * x3 + k[2] / 4.0 * x4;
    return Epot;
}


double computeLinearMomentum(const std::vector<double>& p){
    double P = 0.0;
    for(int i=0; i<p.size(); i++){
        P += p[i];
    }
    return P;
}
double computeCenterOfMass(const std::vector<double>& x){
    double Q = 0.0;
    for(int i=0; i<x.size(); i++){
        Q += x[i];
    }
    return Q;
}



void writeProperties(const std::vector<double>& x, const std::vector<double>& p, const std::vector<double>& k, FILE * outfile){
    double Ekin = computeKineticEnergy(p);
    double Epot = computePotentialEnergy(x, k);
    double meanPsq = computeMeanPsq(p);
    double P = computeLinearMomentum(p);
    double Q = computeCenterOfMass(x);

    fprintf(outfile, "%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", meanPsq, Ekin+Epot, Ekin, Epot, P, Q);
}


void dumpCoordinates(const std::vector<double>& x, const std::vector<double>& p, FILE * coords){
    for (int i=0; i<x.size(); i++){
        fprintf(coords, "%8.4f ", x[i]);
    }
    for (int i=0; i<p.size(); i++){
        fprintf(coords, "%8.4f ", p[i]);
    }
    fprintf(coords, "\n");
}
