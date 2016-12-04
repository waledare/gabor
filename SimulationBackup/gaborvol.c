#include "gaborvol.h"

double complex sigmahat (
        int  H, 
        int  K, 
        double complex * coefs,
        double t) {
    size_t nrow, ncol;
    int i,j,h,k;
    double out;

    nrow = 2*H+1;
    ncol = 2*K+1;

    for (out = 0.00 + 0.00*I, i = 0, h = -H; h <= H; i++, h++)
        for (j = 0, k = -K; k <= K; j++, k++){
            //Because R stores a matrix as an array going down the
            //columns of the matrix first, we need to makeup for it
            //in the indexing of coefs.
            out = out + coefs[j*nrow + i]*ghk(h,k,t);
        }
    return out;
}

void getsigmas (
        int * H, 
        int * K, 
        double * deltaxsq,
        double * realcoef,
        double * imagcoef,
        double * ti,
        double * realsigmas,
        double * imagsigmas,
        size_t * nt){
    int i;
    double complex sigmati;
    double complex * coefs;
    size_t nrow, ncol, ncoef;
    
    if (*H < 0 || *K < 0){
        printf("getsigmaserror: invalid values passed to H or K\n");
        exit(1);
    }

    nrow = 2*H[0]+1;
    ncol = 2*K[0]+1;
    ncoef = nrow*ncol;
 
    coefs = (double complex *) malloc(ncoef*sizeof(double complex));
    if(NULL == coefs){
        printf("getsigmaerror: mem allocation failed for coefs\n");
        exit(1);
    }
    
    //This gets the coefficients.
    gaborcoefficients(H,K,deltaxsq,ti,realcoef,imagcoef,nt);

    tocomplex(realcoef, imagcoef, coefs, ncoef);

    for(sigmati = 0.00 + 0.0*I, i = 0; i < nt[0]; i++){
        sigmati = sigmahat(*H,*K,coefs,ti[i]);
        realsigmas[i] = creal(sigmati);
        imagsigmas[i] = cimag(sigmati);
    }
    free(coefs);
}
