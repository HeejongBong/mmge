#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>

#include "dmatrix.h"
#include "def.h"

extern "C"{
    double MAX_double(double* a, double *b){
        return(a[0] > b[0] ? a[0] : b[0]);
    }
    
    double ABS(double* a){
        return(a[0] > 0 ? a[0] : -a[0]);
    }
    
    void screen_complete(int *AA, double* S_bar, double* lambda1, double* lambda2, int* nn, int* pp, int* LL){
        
        double* S = S_bar;
        int* A = AA;
        double lam1 = lambda1[0];
        double lam2 = lambda2[0];
        int p = pp[0];
        int L = LL[0];
        int* n = nn;
        
        for (int i = 0; i < p; i++) {
            for (int j = 0; j <= i; j++) {
                double* s = (double*) malloc(sizeof(double)*L);
                int index = j*p*L+i;
                for (int l = 0; l < L; l++) {
                    s[l] = S[index+l*p]*n[l];
                }

                std::sort(s, s + L);
                double sum1 = 0;
                double sum2  = 0;
                for (int l = 0; l < L; l++) {
                    sum1 += s[l];
                    sum2 += s[L-l-1];
                    double tmp1 = sum1/(l+1);
                    double tmp2 = sum2/(l+1);
                    tmp1 = ABS(&tmp1);
                    tmp2 = ABS(&tmp2);
                    if (MAX_double(&tmp1,&tmp2) > (lam1 + (L-l)*lam2)) A[i*p+j] = 1;
                }
                if (i == j) {
                    A[i*p+j] = 0;
                }
            }
        }
    }
    
    void test_lapack(){
        double* vec = (double*) malloc(9*sizeof(double));
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                vec[i*3+j] = 1.0+i+j;
            }
        }
        double* ivec1 = (double*) malloc(9*sizeof(double));
        double* ivec2 = (double*) malloc(9*sizeof(double));
        double* ivec3 = (double*) malloc(9*sizeof(double));
        dmat_vset(9,1.0,ivec1);
        dmat_vcopy(9,ivec1,ivec2);
        dmat_waxpby(9,2.0,ivec1,1.0,ivec2,ivec3);
        dmat_elemprod(9,ivec1,ivec2,ivec3);
        dmat_elemdivi(9,ivec1,ivec2,ivec3);
        for (int i = 0; i < 9; i++) {
            printf("%f, ", ivec3[i]);
        }
        printf("\n");
        double norm2 = dmat_norm2(9,vec);
        printf("The l2 norm of vec is %f. \n", norm2);
        
        double* eigvec = (double*) malloc(9*sizeof(double));
        double* eigval = (double*) malloc(3*sizeof(double));
        eigen_decomp(3,vec,eigvec,eigval);
        for (int i = 0; i < 3; i++) {
            printf("%f, ", eigval[i]);
        }
        printf("\n");
        free(vec); free(eigvec); free(eigval);
        
    }
    
    int main(){
        test_lapack();
    }
    
}
