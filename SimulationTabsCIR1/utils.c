#include "utils.h"

int factorial (int n){
    double out;
    int i;
    if ( n < 0){
        printf("factorial error: n is less than 0\n");
        exit(1);
    }
    for (out = 1, i = n; i > 0;i--)
        out = out * i;

    return out ;
}

int combination (int n , int k){
    if ( n < k){
        printf("combination error: n is less than k\n");
        exit(1);
    }
    return factorial(n)/(factorial(k)*factorial(n-k));
}
void tocomplex (
        double * realxs,
        double * imxs,
        double complex * cxs,
        size_t n) {
    size_t i;
    for (i = 0; i < n; i++){
        cxs[i] = realxs[i] + imxs[i]*I;
    }
}

void realtocomplex (double * rs, double complex * cs, size_t n){
    int i;

    for (i = 0; i < n ; i++)
        cs[i] = rs[i] + 0.0*I;
}

double complex csum( double complex * v,  size_t nv){
    int i;
    double complex out;

    for (out = 0.0 + 0.0*I, i = 0; i < nv; i++)
        out = out +  v[i];
    return out;
}

//The second vector always gets the conjugate
double complex cinnerprod( double complex * v1, double complex * v2,  size_t nv){
    int i;
    double complex out;
    double complex * cs;

    if( NULL == (cs = (double complex *)malloc(nv * sizeof(double complex)))){
        printf("cinnerprod error: malloc failed\n");
        exit(1);
    }
    for (i =0; i < nv ; i++)
        cs[i] = v1[i]*conj(v2[i]);
    out = csum( cs, nv);
    free(cs);
    return out;
}
