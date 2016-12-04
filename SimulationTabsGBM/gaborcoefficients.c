#include "gaborcoefficients.h"

void gaborcoefficients (
        int * H, 
        int * K, 
        double * deltaxsq,
        double * ts,
        double * realcoef,
        double * imagcoef,
        int * nts
        ) {
    size_t nrow, ncol, nt;
    int i,j, h, k;
    double complex * gthkti, * cdeltaxsq, ** ccoef;
    
    nrow = 2*H[0]+1;
    ncol = 2*K[0]+1;
    nt = nts[0];
    if (H[0] < 0 || K[0] < 0 ){
        printf("Either H or K is negative\n");
        exit(1);
    }

    if(NULL==(ccoef=(double complex **)malloc(nrow*sizeof(double complex *)))){
        printf("gaborcoefficients error: malloc for ccoef failed\n");
        exit(1);
    }

    for (i = 0; i < nrow ; i++)
        if(NULL == 
                (ccoef[i] =
                 (double complex *)malloc(ncol*sizeof(double complex))) ){
            printf("gaborcoefficients error: malloc for ccoef[i] failed\n");
            exit(1);
       }

    if(NULL == (gthkti = (double complex *)malloc(nt*sizeof(double complex)))){
        printf("gaborcoefficients error: malloc for fhkti failed\n");
        exit(1);
    }

    if(NULL ==(cdeltaxsq=(double complex *)malloc(nt*sizeof(double complex)))){
        printf("gaborcoefficients error: malloc cdeltaxsq failed\n");
        exit(1);
    }
    realtocomplex(deltaxsq, cdeltaxsq, nt);
    for (i = 0, h = -H[0]; h <= H[0]; i++, h++)
        for (j = 0, k = -K[0]; k <= K[0]; j++, k++){
            getgthkti(h,k,ts,gthkti,nt);
            //So R does something weird here. To unfurl the matrix
            //ccoef into an array, it is required to do so going 
            //down the column instead of the row, so we 
            //build the array that way.
            ccoef[i][j] = cinnerprod(cdeltaxsq, gthkti, nt);
            realcoef[j*nrow + i] = creal(ccoef[i][j]);
            imagcoef[j*nrow + i] = cimag(ccoef[i][j]);
        }
    free(gthkti);
    free(cdeltaxsq);
    for (i = 0; i < nrow ; i++)
        free(ccoef[i]);
    free(ccoef);
}

//Get dual frame element values at time ti
void getgthkti (
        int h,
        int k,
        double * ts,
        double complex * gthktis,
        size_t nts){
    int i;

    for(i = 0; i < nts; i++)
        gthktis[i] = gthk(h,k,ts[i]);
}
