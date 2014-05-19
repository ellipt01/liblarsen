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
// Level1
extern double	dasum_  (const int *n, const double *x, const int *incx);
extern void	daxpy_  (const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
extern void	dcopy_  (const int *n, const double *x, const int *incx, double *y, const int *incy);
extern double	ddot_   (const int *n, const double *x, const int *incx, const double *y, const int *incy);
extern double	dnrm2_  (const int *n, const double *x, const int *incx);
extern void	dscal_  (const int *n, const double *alpha, double *x, const int *incx);
extern int		idamax_ (const int *n, const double *x, const int *incx);
// Level2
extern void	dgemv_ (const char *trans, const int *m, const int *n, const double *alpha, const double *a, const int *lda,
			const double *x, const int *incx, const double *beta, double *y, const int *incy);
// Level3
extern void	dgemm_ (const char *transa, const char *transb, const int *m, const int *n, const int *k,
		const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
		const double *beta, double *c, const int *ldc);
#endif

typedef struct s_linsys		linsys;
typedef struct s_penalty		penalty;

/* structure for linear system of regression equations */
struct s_linsys {
	size_t			n;	// number of data
	size_t			p;	// number of variables

	bool			y_centerdized;
	bool			x_centerdized;
	bool			x_normalized;

	double			*y;		// data
	double			*x;		// variables

	double			*meany;
	double			*meanx;
	double			*normx;

	const penalty	*pen;
};

/* penalty term */
struct s_penalty {
	size_t			p1;
	/* scale factor: scale = 1. / (a + b * lambda2) */
	double			a;
	double			b;
	double			*r;
};

/* linsys.c */
linsys		*linsys_alloc (const size_t n, const size_t p, const double *y, const double *x, const double *meany, const double *meanx, const double *normx);
void		linsys_free (linsys *l);
double		*linsys_centering (const size_t size1, const size_t size2, double *x);
double		*linsys_normalizing (const size_t size1, const size_t size2, double *x);

penalty	*penalty_alloc (const size_t p1, const size_t p, const double a, const double b, const double *r);
void		penalty_free (penalty *pen);
void		linsys_set_penalty (linsys *sys, const penalty *pen);

#ifdef __cplusplus
}
#endif

#endif /* LINSYS_H_ */
