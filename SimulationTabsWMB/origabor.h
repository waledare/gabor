#ifndef GABOR_H
#define GABOR_H

#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define a 1
#define b (1.00/3.00)

double g(double x);//Gabor generator
double gt(double x);//Dual Gabor generator
double complex ghk (int h, int k , double x);//Frame elements
double complex gthk (int h, int k , double x);//Dual frame elements
#endif  /*GABOR_H*/
