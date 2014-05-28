/*
 * larsen_private.h
 *
 *  Created on: 2014/05/14
 *      Author: utsugi
 */

#ifndef LARSEN_PRIVATE_H_
#define LARSEN_PRIVATE_H_

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

/* Private macros and entry of constants
 * which are only used internally */

/* cast pointer type to const int * */
#ifndef CINTP
#define CINTP(a)	((const int *) &(a))
#endif

/* min(a, b) */
#ifndef LARSEN_MIN
#define LARSEN_MIN(a, b)	(((a) > (b)) ? (b) : (a))
#endif

/* following constants are set in larsen_linalg.c */
extern const int		ione;	//  1
extern const double	dzero;	//  0.
extern const double	dmone;	// -1.

/* utils.c */
void	larsen_error (const char * function_name, const char *error_msg);

/* linalg.c */
void	larsen_awpy (const larsen *l, double alpha, double *w, double *y);
int		larsen_linalg_cholesky_svx (const size_t size, double *l, const size_t lda, double *b);
int		larsen_linalg_cholesky_insert (const size_t n, double **r, const int index, double *u);
void	larsen_linalg_cholesky_delete (const size_t size, double **r, const int index);

/* active_set.c */
bool	update_activeset (larsen *l);
int		*complementA (larsen *l);

/* stepsize.c */
bool	update_stepsize (larsen *l);

/* equiangular.c */
bool	update_equiangular (larsen *l);

#endif /* LARSEN_PRIVATE_H_ */
