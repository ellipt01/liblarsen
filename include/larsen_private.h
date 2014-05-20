/*
 * larsen_private.h
 *
 *  Created on: 2014/05/20
 *      Author: utsugi
 */

#ifndef LARSEN_PRIVATE_H_
#define LARSEN_PRIVATE_H_

/* Private macros, constants and headers
 * which are only used internally */

#include "linsys_private.h"

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

/* linalg.c */
void		larsen_axapy (larsen *l, double alpha, double *xa, double *y);
double		*larsen_xa_dot_ya (larsen *l, const size_t n, double alpha, const double *x, const double *ya);
double		*larsen_xa_transpose_dot_y (larsen *l, const size_t n, const double alpha, const double *x, const double *y);

/* active_set.c */
int			*complementA (larsen *l);
bool		update_activeset (larsen *l);

/* stepsize.c */
bool		update_stepsize (larsen *l);

/* equiangular.c */
bool		update_equiangular (larsen *l);

#endif /* LARSEN_PRIVATE_H_ */
