#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

#define MAX(a,b) (a > b) ? a : b
int factorial (int n);
int combination (int n, int k);
void tocomplex (
        double * realxs,
        double * imxs,
        double complex * cxs,
        size_t n);
void realtocomplex (
        double * rs, 
        double complex * cs, 
        size_t n);
double complex csum( 
        double complex * v,  
        size_t nv);
double complex cinnerprod( 
        double complex * v1, 
        double complex * v2,  
        size_t nv);
#endif  /*UTILS_H*/
