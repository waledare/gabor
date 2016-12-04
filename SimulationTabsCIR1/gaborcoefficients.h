#ifndef GABORCOEFFICIENTS_H
#define GABORCOEFFICIENTS_H

#include "utils.h"
#include "gabor.h"

void gaborcoefficients (
        int * H, 
        int * K, 
        double * deltaxsq,
        double * ts,
        double * realcoef,
        double * imagcoef,
        int * nts
        ) ;
 
void getgthkti (
        int h,
        int k,
        double * ts,
        double complex * fhktis,
        size_t nts);
 

#endif  /*GABORCOEFFICIENTS_H*/
