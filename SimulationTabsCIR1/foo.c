#include "foo.h"

void foo (
        int* nin, 
        double* realxs, 
        double * imxs) {
    int i;
    double complex * cmplx;
    cmplx = (double complex *) (malloc(nin[0] * sizeof(double complex)));
    if (NULL == cmplx){
        printf("complex array not allocated \n");
        exit(1);
    }
    tocomplex(realxs, imxs, cmplx, (size_t)nin[0]);
    for(i = 0; i < nin[0] ; i++){
        cmplx[i]= (cmplx[i]*cmplx[i]);
        realxs[i] = creal(cmplx[i]);
        imxs[i] = cimag(cmplx[i]);
    }
}

