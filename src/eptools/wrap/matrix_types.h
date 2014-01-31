/* -------------------------------------------------------------------
 * Type and macro definitions related to matrix/vectors and BLAS
 * functions
 * -------------------------------------------------------------------
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#ifndef MATRIX_TYPES_H
#define MATRIX_TYPES_H

/*
 * Types for matrix arguments to support BLAS convention.
 * Matrix is n-by-m, stored column-major (Fortran convention) in 'buff'.
 * Column i starts at 'buff[i*stride]', where stride>=m. 'strcode':
 * - &UPLO(strcode): "L" or "U" (lower or upper triangular)
 * - &DIAG(strcode): "N" or "U" (normal or unit diagonal)
 */
typedef struct {
  double* buff;
  int m,n;
  int stride;
  char strcode[4];
} fst_matrix;

typedef struct {
  double* buff;
  int n;
  int stride;
} fst_vector;

/*
 * Macros picking out specific positions from structure code array
 */
#define UPLO(arr) (arr)[0]

#define DIAG(arr) (arr)[2]

/*
 * The BLAS/LAPACK function xxx is called as xxx_ in Linux, but as
 * xxx in Windows. Uncomment the corresponding definition here.
 */
#ifdef HAVE_WINDOWS
/* Version for Windows */
#define BLASFUNC(NAME) NAME
#else
/* Version for Linux */
#define BLASFUNC(NAME) NAME ## _
#endif

/*
 * Pointer to function types for BLAS functions
 * In 64bit versions of BLAS, 'int' becomes 'ptrdiff_t'.
 */
#ifdef HAVE_BLAS64BIT
#include <stddef.h>
typedef ptrdiff_t blasint_t;
#else
typedef int blasint_t;
#endif

typedef void (* dswap_type) (blasint_t* n,double* x,blasint_t* incx,
			     double* y,blasint_t* incy);

typedef void (* dcopy_type) (blasint_t* n,double* x,blasint_t* incx,double* y,
			     blasint_t* incy);

typedef void (* dscal_type) (blasint_t* n,double* alpha,double* x,
			     blasint_t* incx);

typedef double (* ddot_type) (blasint_t* n,double* a,
			      blasint_t* lda,double* b,
			      blasint_t* ldb);

typedef void (* daxpy_type) (blasint_t *n,double* alpha,
			     double* x,blasint_t* incx,double* y,
			     blasint_t* incy);

typedef void (* dsymv_type) (char* uplo,blasint_t* n,double* alpha,
			     double* a,blasint_t* lda,double* x,
			     blasint_t* incx,double* beta,double* y,
			     blasint_t* incy);

typedef void (* dgemm_type) (char* tra,char* trb,blasint_t* m,
			     blasint_t *n,blasint_t* k,double* alpha,
			     double* a,blasint_t* lda,double* b,
			     blasint_t* ldb,double* beta,double* c,
			     blasint_t* ldc);

typedef void (* dsymm_type) (char* side,char* uplo,blasint_t* m,
			     blasint_t* n,double* alpha,double* a,
			     blasint_t* lda,double* b,
			     blasint_t* ldb,double* beta,double* c,
			     blasint_t* ldc);

typedef void (* dtrsm_type) (char* side,char* uplo,
			     char* trans,char* diag,blasint_t* m,
			     blasint_t* n,double* alpha,double* a,
			     blasint_t* lda,double* b,
			     blasint_t* ldb);

typedef void (* dtrmm_type) (char* side,char* uplo,
			     char* trans,char* diag,blasint_t* m,
			     blasint_t* n,double* alpha,double* a,
			     blasint_t* lda,double* b,
			     blasint_t* ldb);

typedef void (* drotg_type) (double* a,double* b,double* c,double* s);

typedef void (* drot_type) (blasint_t* n,double* x,
			    blasint_t* incx,double* y,
			    blasint_t* incy,double* c,
			    double* s);

typedef void (* dtrsv_type) (char* uplo,char* trans,
			     char* diag,blasint_t* n,double* a,
			     blasint_t* lda,double* x,blasint_t* incx);

#endif
