#include "gabor.h"

//Generator
inline double g (double x) {
    double out;
    if (x < 1.0 && x >= 0.0)
        out = x;
    if (x < 2.0 && x >= 1.0)
        out = 2 - x;
    /*if (x < 1.0 && x >= 0.0)
        out = (x*x)/2.0;
    if (x < 2.0 && x >= 1.0)
        out = (-2*x*x + 6*x - 3)/2.0;
    if (x < 3.0 && x >= 2.0)
        out = ((3-x)*(3-x))/2.0;*/
    return out;
    //return (bspline(SPLINE_ORDER, x));
}

//Dual
double gt (double x) {
    int n;
    double out, a, b ;
    /*if (x < 0.0 && x >= -1.0)
        out = (2.0/3.0)*(x+1);
    if (x < 2.0 && x >= 0.0)
        out = (1.0/3.0)*(2 - x);*/
    
    a = T_PARAM;
    b = E_PARAM;
    /*for (out = 0.0, n = 1; n <SPLINE_ORDER ; n++)
        out = out + g(x+n);*/
    out = a*b*g(x) + 2*a*b*g(x+1);
    //out = a*b*g(x) + 2*a*b*g(x+1)+ 2*a*b*g(x+2);
    return out;
}

//Gabor frame elements
double complex ghk (int h, int k , double x){
    double a,b;
    a = T_PARAM;
    b = E_PARAM;
    return ((1/sqrt(a))*cexp(2*M_PI*I*h*b*x)*g((x - k*a)/a));
}

//Dual Gabor frame elements
double complex gthk (int h, int k , double x){
    double a,b;
    a = T_PARAM;
    b = E_PARAM;
    //return (cexp(2*M_PI*I*h*E_PARAM*x)*gt(x - k*T_PARAM));
    return ((1/sqrt(a))*cexp(2*M_PI*I*h*b*x)*gt((x - k*a)/a));
}

inline double bspline  (int n, double x) {
    int  j ;
    double out, c , p;

    for (out = 0.0, j = 0; j <= n ; j++){
        out = out + pow(-1, j)*combination(n , j)*pow(MAX(0,x-j),n-1);
        p = pow(MAX(0,x-j),n-1);
        c = MAX(0,x-j);
        //if(c > 0 )
          //  printf("%d %d %f %f\n",n,j,p,c);
    }
    out = out *( 1.0/(double)factorial(n - 1));
    return out;
}
