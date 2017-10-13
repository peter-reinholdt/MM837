#include <vector>

void computeForces(const std::vector<double>& x, std::vector<double>& forces, const std::vector<double> k){
    double x1, x2, x3;
    int N = forces.size();
    //zero out old forces
    for(int i=0; i<N; i++){
        forces[i] = 0.0;
    }
    
    //do the middle bits
    for(int i=1; i<N-1; i++){
        x1 = x[i] - x[i+1];
        x2 = x1 * x1;
        x3 = x2 * x1; 
        forces[i] += -(k[0]*x1+k[1]*x2+k[2]*x3);
        x1 = x[i-1] - x[i];
        x2 = x1 * x1;
        x3 = x2 * x1; 
        forces[i] -= -(k[0]*x1+k[1]*x2+k[2]*x3);
    }
    //do the edge forces 
    x1 = x[0] - x[1];
    x2 = x1 * x1;
    x3 = x2 * x1; 
    forces[0] += -(k[0]*x1+k[1]*x2+k[2]*x3);
    x1 = x[N-1] - x[0];
    x2 = x1 * x1;
    x3 = x2 * x1; 
    forces[0] -= -(k[0]*x1+k[1]*x2+k[2]*x3);
    // 
    x1 = x[N-1] - x[0];
    x2 = x1 * x1;
    x3 = x2 * x1; 
    forces[N-1] += -(k[0]*x1+k[1]*x2+k[2]*x3);
    x1 = x[N-2] - x[N-1];
    x2 = x1 * x1;
    x3 = x2 * x1; 
    forces[N-1] -= -(k[0]*x1+k[1]*x2+k[2]*x3);
}
