#ifndef GABORVOL_H
#define GABORVOL_H

#include "gabor.h" 
#include "utils.h" 

double complex sigmahat (
        int  H, 
        int  K, 
        double complex * coefs,
        double t);

void getsigmas (
        int * H, 
        int * K, 
        double * deltaxsq,
        double * realcoef,
        double * imagcoef,
        double * ti,
        double * realsigmas,
        double * imagsigmas,
        size_t * nt);
 
#endif  /*GABORVOL_H*/
