/*
 * linsys.c
 *
 *  Created on: 2014/05/19
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <linsys.h>

#include "linsys_private.h"

const int		ione  =  1;
const double	dzero =  0.;
const double	done  =  1.;
const double	dmone = -1.;

/* print an error message and exit */
void
linsys_error (const char * function_name, const char *error_msg, const char *file, const int line)
{
	fprintf (stderr, "ERROR: %s: %s:%d: %s\n", function_name, file, line, error_msg);
	exit (1);
}

linsys *
linsys_alloc (const double lambda2, const size_t n, const size_t p, const double *y, const double *x)
{
	size_t		np;
	linsys		*lsys;

	if (lambda2 < 0.) linsys_error ("lisys_alloc", "lambda2 must be >= 0.", __FILE__, __LINE__);
	if (!y) linsys_error ("lisys_alloc", "vector *y is empty.", __FILE__, __LINE__);
	if (!x) linsys_error ("lisys_alloc", "matrix *x is empty.", __FILE__, __LINE__);
	if (n <= 0) linsys_error ("lisys_alloc", "n must be > 0.", __FILE__, __LINE__);
	if (p <= 0) linsys_error ("lisys_alloc", "p must be > 0.", __FILE__, __LINE__);

	lsys = (linsys *) malloc (sizeof (linsys));

	lsys->n = n;
	lsys->p = p;
	np = lsys->n * lsys->p;

	lsys->y = (double *) malloc (lsys->n * sizeof (double));
	lsys->x = (double *) malloc (np * sizeof (double));

	dcopy_ (LINSYS_CINTP (n), y, &ione, lsys->y, &ione);
	dcopy_ (LINSYS_CINTP (np), x, &ione, lsys->x, &ione);

	/* By default, data is assumed to be not normalized and not standardized */
	lsys->meany = NULL;
	lsys->meanx = NULL;
	lsys->normx = NULL;
	lsys->ycentered = false;
	lsys->xcentered = false;
	lsys->xnormalized = false;

	lsys->lambda2 = lambda2;

	/* regression type */
	if (lambda2 > dlamch_ ("e")) {	// default: ridge regression
		lsys->regtype = REGULARIZATION_RIDGE;
		lsys->scale2 = 1. / (1. + lambda2);
		lsys->scale = sqrt (lsys->scale2);
	} else {	// lasso
		lsys->regtype = REGULARIZATION_LASSO;
		lsys->scale2 = 1. ;
		lsys->scale = 1.;
	}
	lsys->pen = NULL;

	return lsys;
}

void
linsys_free (linsys *lsys)
{
	if (lsys) {
		if (lsys->y) free (lsys->y);
		if (lsys->x) free (lsys->x);

		if (lsys->meany) free (lsys->meany);
		if (lsys->meanx) free (lsys->meanx);
		if (lsys->normx) free (lsys->normx);

		free (lsys);
	}
	return;
}

/* Centering each column of matrix:
 * x(:, j) -> x(:, j) - mean(x(:, j)) */
static double *
centering (const size_t size1, const size_t size2, double *x)
{
	int		i, j;
	double	*mean = (double *) malloc (size2 * sizeof (double));
	for (j = 0; j < size2; j++) {
		double	*xj = x + LINSYS_INDEX_OF_MATRIX (0, j, size1);
		mean[j] = 0.;
		for (i = 0; i < size1; i++) mean[j] += xj[i];
		mean[j] /= (double) size1;
		for (i = 0; i < size1; i++) xj[i] -= mean[j];
	}
	return mean;
}

/* Normalizing each column of matrix:
 * x(:, j) -> x(:, j) / norm(x(:, j)) */
static double *
normalizing (const size_t size1, const size_t size2, double *x)
{
	int		j;
	double	*nrm = (double *) malloc (size2 * sizeof (double));
	for (j = 0; j < size2; j++) {
		double	alpha;
		double	*xj = x + LINSYS_INDEX_OF_MATRIX (0, j, size1);
		nrm[j] = dnrm2_ (LINSYS_CINTP (size1), xj, &ione);
		alpha = 1. / nrm[j];
		dscal_ (LINSYS_CINTP (size1), &alpha, xj, &ione);
	}
	return nrm;
}

/* centering lsys->y,
 * and set lsys->meany[0] = mean( y ) */
void
linsys_centering_y (linsys *lsys)
{
	lsys->meany = centering (lsys->n, 1, lsys->y);
	lsys->ycentered = true;
	return;
}

/* centering each col of lsys->x,
 * and set lsys->meanx[j] = mean( X(:,j) ) */
void
linsys_centering_x (linsys *lsys)
{
	lsys->meanx = centering (lsys->n, lsys->p, lsys->x);
	lsys->xcentered = true;
	return;
}

/* normalizing each col of lsys->x,
 * and set lsys->normx[j] = mean( X(:,j) ) */
void
linsys_normalizing_x (linsys *lsys)
{
	lsys->normx = normalizing (lsys->n, lsys->p, lsys->x);
	lsys->xnormalized = true;
	return;
}

/* standardizing lsys->x */
void
linsys_standardizing_x (linsys *lsys)
{
	linsys_centering_x (lsys);
	linsys_normalizing_x (lsys);
	return;
}

penalty *
penalty_alloc (const size_t pj, const size_t p, const double *r)
{
	penalty	*pen;
	size_t		t;

	if (!r) linsys_error ("penalty_alloc", "matrix *r is empty.", __FILE__, __LINE__);
	pen = (penalty *) malloc (sizeof (penalty));
	t = pj * p;
	pen->pj = pj;
	pen->p = p;
	pen->r = (double *) malloc (t * sizeof (double));
	dcopy_ (LINSYS_CINTP (t), r, &ione, pen->r, &ione);
	return pen;
}

void
penalty_free (penalty *pen)
{
	if (pen) {
		if (pen->r) free (pen->r);
		free (pen);
	}
	return;
}

/* set lsys->pen = pen, and set lsys->scale2 = 1 / (a + b * lambda2) */
void
linsys_set_penalty (linsys *lsys, const double a, const double b, const penalty *pen)
{
	if (!pen) linsys_error ("linsys_set_penalty", "penalty *pen is empty.", __FILE__, __LINE__);
	if (lsys->p != pen->p)
		linsys_error ("linsys_set_penalty", "penalty *pen->p must be same as linsys *lsys->p.", __FILE__, __LINE__);

	if (lsys->regtype != REGULARIZATION_LASSO) {
		lsys->regtype = REGULARIZATION_USER_DEFINED;
		lsys->scale2 = 1. / (a + b * lsys->lambda2);
		lsys->scale = sqrt (lsys->scale2);
	}
	lsys->pen = pen;
}

/* is lasso? : lambda2 <= eps? */
bool
linsys_is_regtype_lasso (const linsys *lsys)
{
	return (lsys->regtype == REGULARIZATION_LASSO);
}

/* is Ridge? :
 * lambda2 > eps && default lsys->pen == NULL, i.e. default regression type? */
bool
linsys_is_regtype_ridge (const linsys *lsys)
{
	return (lsys->regtype == REGULARIZATION_RIDGE);
}

double
linsys_get_lambda2 (const linsys *lsys)
{
	return lsys->lambda2;
}

double
linsys_get_scale (const linsys *lsys)
{
	return lsys->scale;
}

double
linsys_get_scale2 (const linsys *lsys)
{
	return lsys->scale2;
}

const size_t
linsys_get_n (const linsys *lsys)
{
	return lsys->n;
}

const size_t
linsys_get_p (const linsys *lsys)
{
	return lsys->p;
}

/* return number of rows of penalty matrix: lsys->pen->pj */
const size_t
linsys_get_pj (const linsys *lsys)
{
	if (lsys->pen) return lsys->pen->pj;
	return lsys->p;
}

/* return lsys->x */
const double
*linsys_get_x (const linsys *lsys)
{
	return lsys->x;
}

/* return lsys->y */
const double
*linsys_get_y (const linsys *lsys)
{
	return lsys->y;
}

/* return lsys->pen->r */
const double
*linsys_get_penalty (const linsys *lsys)
{
	if (lsys->pen) return lsys->pen->r;
	return NULL;
}
