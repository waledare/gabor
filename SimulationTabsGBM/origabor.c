#include "gabor.h"

//Generator
double g (double x) {
    double out = 0.0;
    
    if (x < 1.0 && x >= 0.0)
        out = x;
    if (x < 2.0 && x >= 1.0)
        out = 2 - x;
    return out;
}

//Dual
double gt (double x) {
    double out = 0.0;
    
    if (x < 0.0 && x >= -1.0)
        out = (2.0/3.0)*(x+1);
    if (x < 2.0 && x >= 0.0)
        out = (1.0/3.0)*(2 - x);
    return out;
}

//Gabor frame elements
double complex ghk (int h, int k , double x){
    return (cexp(2*M_PI*I*h*b*x)*g(x - k));
}

//Dual Gabor frame elements
double complex gthk (int h, int k , double x){
    return (cexp(2*M_PI*I*h*b*x)*gt(x - k));
}


