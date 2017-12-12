#include <cmath>
#include <iostream>


inline void to_interval(double& x){
    x = fmod(x, M_PI);
}


int main(){
    double x;

    for (int i=-10; i<10; i++){
        x = i;
        std::cout << x << " ";
        to_interval(x);
        std::cout << x << std::endl;
    }
}
