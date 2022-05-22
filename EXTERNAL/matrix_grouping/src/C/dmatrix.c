/** \file   dmatrix.c
 *  \brief  Source file for matrix and vector manipulation functions.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#ifndef Rpackage
#include "blas.h"
#include "lapack.h"
#else
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#endif

#include "def.h"
#include "dmatrix.h"

/****************************************************************************/
/*                                                                          */
/*                           VECTOR OPERATIONS                              */
/*                                                                          */
/****************************************************************************/

/** \brief \f$ \|x\|_1 \f$
 *
 *  Returns 1-norm of a vector x.
 *
 *  @param  n   length of a vector x.
 *  @param  x   pointer to a vector x.
 *  @return     result.
 */
double dmat_norm1(const int n, const double *x)
{
    return F77_CALL(dasum)(&n, x, &ione);
}


/** \brief \f$ \|x\|_2 \f$
 *
 *  Returns 2-norm of a vector x.
 *
 *  @param  n   length of a vector x.
 *  @param  x   pointer to a vector x.
 *  @return     result.
 */
double dmat_norm2(const int n, const double *x)
{
    return F77_CALL(dnrm2)(&n, x, &ione);
}


/** \brief \f$ \|x\|_{\infty} \f$
 *
 *  Returns infinity-norm of a vector x
 *
 *  @param  n   length of a vector x
 *  @param  x   pointer to a vector x
 *  @return     result.
 */
double dmat_norminf(const int n, const double *x)
{
    return fabs(x[F77_CALL(idamax)(&n, x, &ione)-1]);
}


/** \brief \f$ x^Ty \f$
 *
 *  Returns dot product of a vector x and y.
 *
 *  @param  n   length of a vector x.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a vector y.
 *  @return     result.
 */
double dmat_dot(const int n, const double *x, const double *y)
{
    return F77_CALL(ddot)(&n, x, &ione, y, &ione);
}


/** \brief \f$ \mbox{dst}_i \leftarrow \mbox{val} \f$
 *
 *  Sets all the elements of a vector with a constant value.
 *
 *  @param  n   length of a vector.
 *  @param  val constant value to set.
 *  @param  dst pointer to a vector.
 */
void dmat_vset(int n, const double val, double *dst)
{
    while (n-- != 0)
        *dst++ = val;
}

void dmat_iset(int n, const int val, int *dst)
{
    while (n-- != 0)
        *dst++ = val;
}


/** \brief \f$ \mbox{dst} \leftarrow \mbox{src} \f$
 *
 *  Copies a vector.
 *
 *  @param  n   length of vectors.
 *  @param  src pointer to a source vector.
 *  @param  dst pointer to a destination vector.
 */
void dmat_vcopy(const int n, const double *src, double *dst)
{
    F77_CALL(dcopy)(&n, src, &ione, dst, &ione);
}

void dmat_icopy(const int n, const int *src, int *dst)
{
    memcpy(dst, src, sizeof(int)*n);
}


/** \brief \f$ y_i = \exp(x_i) \f$
 *
 *  Computes elementwise exp() of a vector.
 *
 *  @param  n   length of vectors.
 *  @param  x   pointer to a source vector.
 *  @param  y   pointer to a destination vector.
 */
void dmat_yexpx(const int n, const double *x, double *y)
{
#ifdef MKL 
    vdExp(n, x, y);
#else
    #ifdef ACML
        vrda_exp(n, x, y);
    #else
        const double *xi;
        double *yi;
        xi = x+n-1;
        yi = y+n-1;
        do {
            *yi-- = exp(*xi--);
        } while (xi >= x);
    #endif
#endif
}
    // y = A*x //
void dmat_yAx(int m, int n, const double *A, const double *x, double *y)
{
    F77_CALL(dgemv)("T",&n,&m,&done,A,&n,x,&ione,&dzero,y,&ione);
}


/** \brief \f$ y_i = x_i^{1/2} \f$
 *
 *  Computes elementwise sqrt() of a vector.
 *
 *  @param  n   length of vectors.
 *  @param  x   pointer to a source vector.
 *  @param  y   pointer to a destination vector.
 */
void dmat_ysqrtx(const int n, const double *x, double *y)
{
#ifdef MKL 
    vdSqrt(n, x, y);
#else
    const double *xi;
    double *yi;
    xi = x+n-1;
    yi = y+n-1;
    do {
        *yi-- = sqrt(*xi--);
    } while (xi >= x);
#endif
}


/** \brief \f$ y_i = 1/x_i \f$
 *
 *  Computes elementwise inv() of a vector.
 *
 *  @param  n   length of vectors.
 *  @param  x   pointer to a source vector.
 *  @param  y   pointer to a destination vector.
 */
void dmat_yinvx(const int n, const double *x, double *y)
{
#ifdef MKL 
    vdInv(n, x, y);
#else
    const double *xi;
    double *yi;
    xi = x+n-1;
    yi = y+n-1;
    do {
        *yi-- = 1/(*xi--);
    } while (xi >= x);
#endif
}


/** \brief \f$ w = \alpha*x + \beta*y \f$
 *
 *  Computes weighted vector sum.
 *  - w = -x
 *  - w = alpha*x
 *  - w = -x + y
 *  - w = alpha*x + y
 *  - w = -x - y
 *  - w = alpha*x - y
 *  - w = alpha*x + beta*y
 *
 *  @param  n       length of vectors.
 *  @param  alpha   constant
 *  @param  x       pointer to a vector.
 *  @param  beta    constant
 *  @param  y       pointer to a vector.
 *  @param  w       pointer to a result vector.
 */
void dmat_waxpby(int n, double alpha, const double *x, double beta,
                 const double *y, double *w)

{
#if 1
    if (w != x && w != y)
    {
        dmat_vset(n, 0, w);
        F77_CALL(daxpy)(&n, &alpha, x, &ione, w, &ione);
        F77_CALL(daxpy)(&n, &beta , y, &ione, w, &ione);
    }
    else if (w == x && w == y)
    {
        double tmp;
        tmp = alpha+beta;
        F77_CALL(dscal)(&n, &tmp, w, &ione);
    }
    else if (w == x /*&& w != y */)
    {
        if (alpha != 1.0) F77_CALL(dscal)(&n, &alpha, w, &ione);
        if (beta  != 0.0) F77_CALL(daxpy)(&n, &beta , y, &ione, w, &ione);
    }
    else /* if (w == y && w != x ) */
    {
        if (beta  != 1.0) F77_CALL(dscal)(&n, &beta , w, &ione);
        if (alpha != 0.0) F77_CALL(daxpy)(&n, &alpha, x, &ione, w, &ione);
    }
#else
    int i;

    if (beta == 0.0)
    {
        if (alpha == -1.0)
        {
            for (i = 0; i < n; i++)
                w[i] = -x[i];
        }
        else
        {
            for (i = 0; i < n; i++)
                w[i] = alpha*x[i];
        }
    }
    else if (beta == 1.0)
    {
        if (alpha == -1.0)
        {
            for (i = 0; i < n; i++)
                w[i] = -x[i] + y[i];
        }
        else
        {
            for (i = 0; i < n; i++)
                w[i] = alpha*x[i] + y[i];
        }
    }
    else if (beta == -1.0)
    {
        if (alpha == -1.0)
        {
            for (i = 0; i < n; i++)
                w[i] = -x[i] - y[i];
        }
        else
        {
            for (i = 0; i < n; i++)
                w[i] = alpha*x[i] - y[i];
        }
    }
    else
    {
        for (i = 0; i < n; i++)
            w[i] = alpha*x[i] + beta*y[i];
    }
#endif
}


/**  \brief \f$ z_i = x_i*y_i \f$
 *
 *  Computes elementwise product of vectors.
 *  
 *  NOTE: x = x.*y is not allowed, i.e., w should not be the same with x.
 *
 *  @param  n   length of a vector x.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a vector y.
 *  @param  z   pointer to a result vector z.
 */
void dmat_elemprod(const int n, const double *x, const double *y, double *z)
{
    if (y != z) F77_CALL(dcopy)(&n, y, &ione, z, &ione);
    F77_CALL(dtbmv)("U", "N", "N", &n, &izero, x, &ione, z, &ione);
}

/** \brief \f$ z_i = x_i/y_i \f$
 *
 *  Computes elementwise division of vectors.
 *
 *  NOTE: y = x./y is not allowed, i.e., w should not be the same with y.
 *
 *  @param  n   length of a vector x.
 *  @param  x   pointer to a vector x.
 *  @param  y   pointer to a vector y.
 *  @param  z   pointer to a result vector z.
 */
void dmat_elemdivi(const int n, const double *x, const double *y, double *z)
{
#ifdef MKL
   vdDiv(n, x, y, z);
#else 
    if (x != z) F77_CALL(dcopy)(&n, x, &ione, z, &ione);
    F77_CALL(dtbsv)("U", "N", "N", &n, &izero, y, &ione, z, &ione);
#endif
}


/****************************************************************************/
/*                                                                          */
/*                           MATRIX OPERATIONS                              */
/*                                                                          */
/****************************************************************************/

/*
 This function calculates the eigenvalues and eigenvectors of
 the n*n symmetric matrix X.
 The matrices have to be in Fortran vector format.
 The eigenvectors will be put columnwise in the n*n matrix eigvec,
 where the corresponding eigenvalues will be put in the vector
 eigval (length n of course). Only the lower triangle of the matrix
 X is used. The content of X is not changed.
 
 This function first queries the Lapack routines for optimal workspace
 sizes. These memoryblocks are then allocated and the decomposition is
 calculated using the Lapack function "dsyevr". The allocated memory
 is then freed.
 */


void dmat_yATx(int m, int n, const double *A, const double *x, double *y)
{
    F77_CALL(dgemv)("N",&n,&m,&done,A,&n,x,&ione,&dzero,y,&ione);
}


/** \brief \f$ y = x'Ax \f$
 *
 *  Computes quadratic forms 
 *
 *  @param  A   pointer to a matrix.
 *  @param  x   pointer to a vector x.
 */
double dmat_xAx(int n, const double *A, const double *x)
{
    double y[n];
    F77_CALL(dgemv)("T",&n,&n,&done,A,&n,x,&ione,&dzero,y,&ione);
    return dmat_dot(n, x, y);
}

/** \brief \f$ C = A^TB \f$
 *
 *  Computes dense matrix-transpose-matrix product.
 *
 *  @param  A       pointer to m by n1 matrix.
 *  @param  B       pointer to m by n2 matrix.
 *  @param  C       pointer to n1 by n2 matrix.
 */
void dmat_C_ATB(int m, int n1, int n2, double* A, double* B, double* C)
{
    F77_CALL(dgemm)("N","T",&n2,&n1,&m,&done,B,&n2,A,&n1,&dzero,C,&n2);
}

/** \brief \f$ B = A^TA \f$
 *
 *  Computes dense matrix-transpose-matrix product.
 *
 *  @param  A       pointer to m by n matrix.
 *  @param  B       pointer to a result matrix.
 */
void dmat_B_ATA(int m, int n, double *A, double *B)
{
    F77_CALL(dsyrk)("L","N",&n,&m,&done,A,&n,&dzero,B,&n);
}

void eigen_decomp(int n, double* X, double *eigvec, double *eigval) {
    
    double *WORK;
    double abstol, WORKopt, vl, vu;
    int *IWORK;
    int numeig, sizeWORK, sizeIWORK, IWORKopt, il, iu,info;
    vl = 0.0;
    vu = 0.0;
    il = 0;
    iu = 0;
    /*  The support of the eigenvectors. We will not use this but the routine needs it  */
    int ISUPPZ[2*n];
    abstol = -1.0; // default tolerance
    
    /*  Query the Lapack routine for optimal sizes for workspace arrays  */
    sizeWORK = -1;
    sizeIWORK = -1;
    F77_CALL(dsyevr)("V", "A", "L", &n, X, &n, &vl,&vu,&il,&iu, &abstol, &numeig, eigval, eigvec, &n, ISUPPZ, &WORKopt, &sizeWORK, &IWORKopt, &sizeIWORK,&info);
    sizeWORK = (int)WORKopt;
    sizeIWORK = IWORKopt;
    WORK = (double*)malloc (sizeWORK*sizeof(double));
    IWORK = (int*)malloc (sizeIWORK*sizeof(int));
    /*  Now calculate the eigenvalues and vectors using optimal workspaces  */
    F77_CALL(dsyevr)("V", "A", "L", &n, X, &n, &vl,&vu,&il,&iu, &abstol, &numeig, eigval, eigvec, &n, ISUPPZ, WORK, &sizeWORK, IWORK, &sizeIWORK,&info);
    /*  Cleanup  */
    free((void*)(WORK)); free((void*)IWORK);
}



