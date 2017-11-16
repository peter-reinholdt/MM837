#include <vector>
#include <iostream>
#include <fstream>
#include "forces.h"

double simpsons_rule(const std::vector<double>& f, const double& delta) { 
	double ret=0.0;
	for (int i = 0 ; i<f.size()-3; i++){
		ret += f[i] + 4.0*f[i+1] + f[i+2];
    }
	ret *= delta/3.0;
	return ret;
}


template<class Potential>
class wavefunction{
    private:
        std::vector<double> psi, psiprime;
        std::vector<double> psiprimeprime;
        std::vector<double> x, x2, ones, V;
        const Potential& pot;
        double xmin, xmax, epsilon, E;
        int nsteps;
    public:
        //we assume normalized wavefunction
        wavefunction(const std::vector<double>& psi_in, 
                     const std::vector<double>& psiprime_in,
                     Potential pot_in,
                     double xmin_in,
                     double xmax_in,
                     int nsteps_in,
                     double E_in):
                     psi(psi_in), psiprime(psiprime_in), pot(pot_in), xmax(xmax_in), xmin(xmin_in), nsteps(nsteps_in), E(E_in) {
                        epsilon= (xmax-xmin)/double(nsteps-1);
                        x.resize(psi.size());
                        x2.resize(psi.size());
                        ones.resize(psi.size());
                        V.resize(psi.size());
                        psiprimeprime.resize(psi.size());
                        double pos = xmin;
                        for (int i=0; i<psi.size(); i++){
                            x[i] = pos;
                            x2[i] = pos * pos;
                            ones[i] = 1.0;
                            V[i] = pot_in.get_pot(pos);
                            pos += epsilon;
                        }
                        //second derivative of psi or first derivative of psiprime.
                        //Do central difference in the middle bits, forward/backward on ends
                        psiprimeprime[0] = (psiprime[1] - psiprime[0])/(epsilon);
                        psiprimeprime[psi.size()-1] = (psiprime[psi.size()-1] - psiprime[psi.size()-2])/(epsilon);
                        for (int i=1; i<psi.size()-1; i++){
                            psiprimeprime[i] = (psiprime[i+1] - psiprime[i-1]) / (2*epsilon);
                        }
                     }


        double expectation_value(const std::vector<double>& left, const std::vector<double>& Operator, const std::vector<double>& right){
            std::vector<double> expectation;
            expectation.resize(left.size());
            for (int i=0; i<left.size(); i++){
                expectation[i] = left[i]*Operator[i]*right[i];
            }
            return simpsons_rule(expectation, epsilon);
        }
    
        //<0|0>
        double norm(){
            return expectation_value(psi, ones, psi);
        }

        //<0|x|0>
        double x_expectation(){
            return expectation_value(psi, x, psi);
        }

        //<0|xx|0>
        double x2_expectation(){
            return expectation_value(psi, x2, psi);
        }

        //<0|p|0>
        //no factor of -i*hbar
        double p_expectation(){
            return expectation_value(psi, ones, psiprime);
        }

        //<0|p2|0>=<0|-1.0*d2/dx2|0>
        //missing factor -1.0*hbar**2
        double p2_expectation(){
            return expectation_value(psi, ones, psiprimeprime);
        }

        //<0|V|0>
        double V_expectation(){
            return expectation_value(psi, V, psi);
        }
        
        // E[x**2] - E[x]**2
        double sigma_x(){
            return x2_expectation() - x_expectation() * x_expectation();
        }

        // E[p**2] - E[p]**2
        double sigma_p(){
            //factors due to not being included earlier
            //hbar = 1                                        (-i) * (-i) = i*i = -1
            return -1.0*p2_expectation() - (-1.0) * p_expectation() * p_expectation(); 
        }
    void write_data(std::string filename){
        std::ofstream fout(filename.c_str());
        for (int i=0; i<psi.size(); i++){
            fout << x[i] << "," << psi[i] << "," << V[i] << "," << E <<"\n";
        } 
    }
};
