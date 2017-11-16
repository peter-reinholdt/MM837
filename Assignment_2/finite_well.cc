/****************************************************************************
 * finite_well.cc
 * obtain energies for finite square well
 * **************************************************************************/

#include<iostream>
#include<iomanip>
#include<sstream>
#include<fstream>
#include "forces.h"
#include "integrators.h"
#include "wavefunction.h"
#include "stdio.h"


using namespace std;


int main(int argc, char** argv) {

	if (argc != 5) {
		cerr << "usage: " << argv[0] <<" <e> <niter> <V0> <a>" << endl;
		return 1;
	}
	double e0; stringstream(argv[1]) >> e0;  //Starting energy. 
	int niter; stringstream(argv[2]) >> niter;  //The number of steps. 
    double V0; stringstream(argv[3]) >> V0;
    double a; stringstream(argv[4]) >> a;

	//Fix the values of phi,phi' at left and right integration endpoints. 
	//These endpoints are [xmin,xmax]
	
    double xmin=-10; 
	double xmax=10;

	double q0_l=exp(-0.5*xmin);
	double q0_r=exp(-0.5*xmax);
	
	double p0_l=0.0;
	double p0_r=0.0;
    
    bool do_nr = true;

	double tol = 1e-6;
	double ecur=e0;
    //double V0 = 16.03074550;
    //double a  = 1.0;
	double de=0.01;
	vector<double> qn_left(1,q0_l), pn_left(1,p0_l);
	vector<double> qn_right(1,q0_r), pn_right(1,p0_r);
	
	int nmax = 1e3;


	if ((niter%2)==0)
		niter++;
	
	double eps=(xmax-xmin)/double(niter-1);

	//The matching point. 
	int imatch = niter/2+10;

	//number of integration steps required from the left
	//and right.
	int niter_l = imatch+1;
	int niter_r = niter-imatch;

	//The results for the wavefunction 
	vector<double> phi(niter, 0.0); 
	vector<double> phi_prime(niter, 0.0); 

	string phi_file("phi.dat");
	ofstream fout(phi_file.c_str());

	double fold;
    double fprime;

	//--------------------------------------------------------------------------

    int it =0; 
	try {

		for (it=0 ; it< nmax; it++) {  
			FiniteSquareWell sw(ecur, V0, a); 
	
			RK4Integrator<FiniteSquareWell> integrator_left(sw, eps);
			integrator_left.setInitialConditions(qn_left,pn_left); 
			
			RK4Integrator<FiniteSquareWell> integrator_right(sw, -eps);
			integrator_right.setInitialConditions(qn_right,pn_right); 

			phi[0] = q0_l;
			phi_prime[0] = p0_l;

			phi[niter-1] = q0_r;
			phi_prime[niter-1] = p0_r; 
			for (int n=1; n<niter_l; n++) { 
				integrator_left(qn_left,pn_left,xmin+(n-1)*eps);
				
				phi[n] = qn_left[0];
				phi_prime[n] = pn_left[0];
			}
			
			for (int n=1; n<niter_r; n++) { 
				integrator_right(qn_right,pn_right,xmax-(n-1)*eps);
				
				phi[niter-n-1] = qn_right[0];
				phi_prime[niter-n-1] = pn_right[0];
			}
			
			//Normalize the left integration results
			for (int n=0; n<niter_l-1; n++) { 
				phi[n] *= qn_right[0]/qn_left[0];
				phi_prime[n] *= qn_right[0]/qn_left[0];
			}

			pn_left[0] *= qn_right[0]/qn_left[0]; 
			double f = (pn_right[0] - pn_left[0])/(pn_right[0] + pn_left[0]);
			
			cout <<setprecision(10)<< "ecurr = "<<ecur <<", pn_right = " << pn_right[0] << ", pn_left = "<<
				pn_left[0] << ", f = "<< f<< endl;

			if (abs(f) < tol) 
				break;
            
            if(do_nr){
                //Newton-Rhapson
                if (it>0){
                    if (fabs(fold) < fabs(f)){
                        //backtrack
                        ecur -= de;
                        fprime = (f-fold)/(de);
                        de = 0.1* (-f/fprime);
                    }else{
                    fprime = (f-fold)/(de);
                    de =  -f/fprime;
                    }
                }
            }else{
                if (it>0){
                    if (abs(fold)<abs(f)){
                        de = -0.5*de; 
                    }
                }
            }
            

			ecur += de; 
			fold = f;
			qn_left[0]=q0_l; pn_left[0]=p0_l;
			qn_right[0]=q0_r; pn_right[0]=p0_r;

		}

	} catch (const invalid_argument& e) {
		cerr << "Invalid argument: "<< e.what() << endl;
		return 1;
	}
	catch (const logic_error& e) {
		cerr << "Logic Error: "<< e.what() << endl;
		return 2;
	}

	if (it==nmax)
		cout << "WARNING: max iterations reached" << endl;

	cout << setprecision(9) << "it = "<<it<<", e = " <<ecur<<", de = " << de<< ", f = " <<fold<<endl;
	
	cout << "Normalizing..." << endl;

	vector<double> norm_arr; 
	for (int i = 0; i<phi.size(); i++) 
		norm_arr.push_back(phi[i]*phi[i]);

	double norm = simpsons_rule(norm_arr, eps); 
	cout << "norm = " << norm << endl;	





	fout << setprecision(6); 
	for (int i = 0 ; i<phi.size(); i++) {
		phi[i] /= sqrt(norm);
		phi_prime[i] /= sqrt(norm);

		fout << (double)(i)*eps+xmin << "     " << phi[i] << endl;
	}


	vector<double> norm_arr2; 
	for (int i = 0; i<phi.size(); i++) 
		norm_arr2.push_back(phi[i]*phi[i]);

	norm = simpsons_rule(norm_arr2, eps); 
	cout << "norm = " << norm << endl;	
    
    FiniteSquareWell sw(ecur, V0, a); 
    wavefunction<FiniteSquareWell> wfn(phi, phi_prime, sw, xmin, xmax, niter, ecur);
    cout << "Norm = " << wfn.norm() << "\n";
    cout << "<x> = " << wfn.x_expectation() << "\n";
    cout << "<x2> = " << wfn.x2_expectation() << "\n";
    cout << "<p> = " << wfn.p_expectation() << "\n";
    cout << "<p2> = " << wfn.p2_expectation() << "\n";
    cout << "sigma(p)**2 = " << wfn.sigma_p()<< "\n";
    cout << "sigma(x)**2 = " << wfn.sigma_x()<< "\n";
    cout << "sigma(x)*sigma(p) = " << sqrt(wfn.sigma_p() * wfn.sigma_x()) << "\n";
    cout << "<T> = " << -1.0 * wfn.p2_expectation() << "\n";
    cout << "<V> = " << wfn.V_expectation() << "\n";
    char buffer [256];
    sprintf(buffer, "fw_E_%e_V_%e_a_%e.dat", ecur, V0, a);
    wfn.write_data(buffer);
	return 0;

}
