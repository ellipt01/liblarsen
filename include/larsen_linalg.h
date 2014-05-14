/*
 * larsen_linalg.h
 *
 *  Wrapper of lapack and qrupdate
 *
 *  Created on: 2014/04/08
 *      Author: utsugi
 */

#ifndef LARSEN_LINALG_H_
#define LARSEN_LINALG_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef LARSEN_INDEX_OF_MATRIX
#define LARSEN_INDEX_OF_MATRIX(i, j, lda) ((i) + (j) * (lda))
#endif

/* blas */
#ifdef HAVE_BLAS_H
#include <blas.h>
#else
// level 1
extern double	dasum_  (const int *n, const double *x, const int *incx);
extern void	daxpy_  (const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
extern void	dcopy_  (const int *n, const double *x, const int *incx, double *y, const int *incy);
extern double	ddot_   (const int *n, const double *x, const int *incx, const double *y, const int *incy);
extern double	dnrm2_  (const int *n, const double *x, const int *incx);
extern void	dscal_  (const int *n, const double *alpha, double *x, const int *incx);
extern int		idamax_ (const int *n, const double *x, const int *incx);
// level 2
extern void	dgemv_ (const char *trans, const int *m, const int *n, const double *alpha, const double *a, const int *lda,
			const double *x, const int *incx, const double *beta, double *y, const int *incy);
#endif

/* lapack */
#ifdef HAVE_LAPACK_H
#include <lapack.h>
#else
extern double	dlamch_ (const char *cmach);
extern void	dpotrs_ (const char *uplo, const int *n, const int *nrhs, const double *a, const int *lda, double *b, const int *ldb, int *info);
#endif

/* qrupdate: choesky linsert/delete */
#ifdef HAVE_QRUPDATE_H
#include <qrupdate.h>
#else
extern void	dchinx_ (const int *n, double *R, const int *ldr, const int *j, const double *u, double *w, int *info);
extern void	dchdex_ (const int *n, double *R, const int *ldr, const int *j, double *w);
#endif

int		larsen_linalg_cholesky_svx (const size_t size, double *l, const size_t lda, double *b);
int		larsen_linalg_cholesky_insert (const size_t n, double **r, const int index, double *u);
void	larsen_linalg_cholesky_delete (const size_t size, double **r, const int index);

#ifdef __cplusplus
}
#endif

#endif /* LARSEN_LINALG_H_ */
