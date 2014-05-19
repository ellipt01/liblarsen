/*
 * linsys.h
 *
 *  Wrapper of lapack and qrupdate
 *
 *  Created on: 2014/04/08
 *      Author: utsugi
 */

#ifndef LINSYS_H_
#define LINSYS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdbool.h>

#ifndef LINSYS_INDEX_OF_MATRIX
#define LINSYS_INDEX_OF_MATRIX(i, j, lda) ((i) + (j) * (lda))
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

typedef struct s_linsys	linsys;

/* structure for linear system of regression equations */
struct s_linsys {
	size_t		n;	// number of data
	size_t		p;	// number of variables

	bool		y_centerdized;
	bool		x_centerdized;
	bool		x_normalized;

	double		*y;		// data
	double		*x;		// variables

	double		*meany;
	double		*meanx;
	double		*normx;

};

/* linsys.c */
linsys		*linsys_alloc (const size_t n, const size_t p, const double *y, const double *x, const double *meany, const double *meanx, const double *normx);
void		linsys_free (linsys *l);
double		*linsys_centering (const size_t size1, const size_t size2, double *x);
double		*linsys_normalizing (const size_t size1, const size_t size2, double *x);

#ifdef __cplusplus
}
#endif

#endif /* LINSYS_H_ */
