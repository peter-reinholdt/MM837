/*****************************************************************************
 * forces.h 
 *
 * Some functors which can be used with the Integrator classes in integrators.h
 *****************************************************************************/

#ifndef __examples_h__
#define __examples_h__

#include <cmath>
#include <vector>
#include <iostream>


class InfiniteSquareWell { 
	double e;
	public: 
		InfiniteSquareWell( const double& e_in) : e(e_in) {} 
		int getNDof() const { 
			return 1;
		}
		double operator() (int i, const std::vector<double>& q, 
				const double& x) const { 
			return -e*q[0]; 
		}


        double get_pot(double x){
            return 0.0;
        }
};


class FiniteSquareWell{ 
	double e, V0, a;
	public: 
		FiniteSquareWell(const double& e_in, const double& V0_in, const double& a_in) : e(e_in), V0(V0_in), a(a_in) {} 
		int getNDof() const { 
			return 1;
		}
		double operator() (int i, const std::vector<double>& q, const double& x) const { 
            double v;
            if (fabs(x) <= a){
                v = -V0;
            }else{
                v = 0.0;
            }
			return (v-e)*q[0]; 
		}


        double get_pot(double x){
            double v;
            if (fabs(x) <= a){
                v = -V0;
            }else{
                v = 0.0;
            }
            return v;
        }
};


class LinearWell{ 
	double e, V0, a;
	public: 
		LinearWell(const double& e_in, const double& V0_in, const double& a_in) : e(e_in), V0(V0_in), a(a_in) {} 
		int getNDof() const { 
			return 1;
		}
		double operator() (int i, const std::vector<double>& q, const double& x) const { 
            double v;
            if (fabs(x) <= a){
                v = V0/a * (fabs(x) - a);
            }else{
                v = 0.0;
            }
			return (v-e)*q[0]; 
		}


        double get_pot(double x){
            double v;
            if (fabs(x) <= a){
                v = V0/a * (fabs(x) - a);
            }else{
                v = 0.0;
            }
            return v;
        }
};



class HarmonicOscillator { 
	double e, k;
	public: 
		HarmonicOscillator( const double& e_in, const double& k_in): e(e_in),
	    k(k_in) {} 
		int getNDof() const { 
			return 1;
		}
		double operator() (int i, const std::vector<double>& q, 
				const double& x) const {
			double v = k*x*x; 
			return (v-e)*q[0]; 
		}


        double get_pot(double x){
            double v = k * x * x;
            return v;
        }
};
#endif 
