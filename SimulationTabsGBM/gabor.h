#ifndef GABOR_H
#define GABOR_H

#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <tgmath.h>
#include "spline.h"

#define T_PARAM 1.0/8.0
#define SPLINE_ORDER 2
#define E_PARAM (1.00/(2*SPLINE_ORDER -1 + 2))

inline double g(double x);//Gabor generator
double gt(double x);//Dual Gabor generator
double complex ghk (int h, int k , double x);//Frame elements
double complex gthk (int h, int k , double x);//Dual frame elements
inline double bspline  (int n, double x) ;
#endif  /*GABOR_H*/
