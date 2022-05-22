#ifndef DMATRIX_H
#define DMATRIX_H
/** \file   dmatrix.h
 *  \brief  Header file for matrix and vector manipulation functions.
 */

#ifdef __cplusplus
extern "C" {
#endif


double dmat_norm1(const int n, const double *x);
double dmat_norm2(const int n, const double *x);
double dmat_norminf(const int n, const double *x);

double dmat_dot(const int n, const double *x, const double *y);

void dmat_vset(int n, const double val, double *dst);
void dmat_iset(int n, const int val, int *dst);
void dmat_vcopy(const int n, const double *src, double *dst);
void dmat_icopy(const int n, const int *src, int *dst);
void dmat_yexpx(const int n, const double *x, double *y);
void dmat_ysqrtx(const int n, const double *x, double *y);
void dmat_yinvx(const int n, const double *x, double *y);
void dmat_yAx(int m, int n, const double *A, const double *x,double *y);
void dmat_waxpby(int n, double alpha, const double *x, double beta,
                 const double *y, double *w);
void dmat_C_ATB(int m, int n1, int n2, double* A, double* B, double* C);

void dmat_elemprod(const int n, const double *x, const double *y, double *z);
void dmat_elemdivi(const int n, const double *x, const double *y, double *z);
void dmat_B_ATA(int m, int n, double *A, double *B);
void dmat_yATx(int m, int n, const double *A, const double *x, double *y);


/* print for debugging */
void dmat_vprint(const int n, const double *v);
double dmat_xAx(int n, const double *A, const double *x);
void eigen_decomp(int n, double* X, double *eigvec, double *eigval);


#ifdef __cplusplus
}
#endif

#endif /* DMATRIX_H */
