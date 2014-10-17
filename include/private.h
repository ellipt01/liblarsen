/*
 * private.h
 *
 *  Created on: 2014/05/14
 *      Author: utsugi
 */

#ifndef PRIVATE_H
#define PRIVATE_H

/* Private macros, constants and headers
 * which are only used internally */

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
extern void	dsymv_ (const char *uplo, const int *n, const double *alpha, const double *a, const int *lda,
		const double *x, const int *incx, const double *beta, double *y, const int *incy);
// Level3
extern void	dgemm_ (const char *transa, const char *transb, const int *m, const int *n, const int *k,
		const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
		const double *beta, double *c, const int *ldc);
#endif

/* lapack */
#ifdef HAVE_LAPACK_H
#include <lapack.h>
#else
extern void	dpotrs_ (const char *uplo, const int *n, const int *nrhs, const double *a, const int *lda, double *b, const int *ldb, int *info);
#endif

/* qrupdate: choesky linsert/delete */
#ifdef HAVE_QRUPDATE_H
#include <qrupdate.h>
#else
extern void	dchinx_ (const int *n, double *R, const int *ldr, const int *j, const double *u, double *w, int *info);
extern void	dchdex_ (const int *n, double *R, const int *ldr, const int *j, double *w);
#endif

/* following constants are set in private.c */
extern const int		ione;	//  1
extern const double	dzero;	//  0.
extern const double	done;	//  1.
extern const double	dmone;	// -1.

/* print error message and terminate program */
void	error_and_exit (const char * function_name, const char *error_msg, const char *file, const int line);
/* print warning message */
void	printf_warning (const char * function_name, const char *error_msg, const char *file, const int line);

/* Private macros and entry of constants
 * which are only used internally */

/* min(a, b) */
#ifndef LARSEN_MIN
#define LARSEN_MIN(a, b)	(((a) > (b)) ? (b) : (a))
#endif

#define CINTP(a) ((const int *) &(a))

#ifdef DEBUG
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
#endif

#endif /* PRIVATE_H */
