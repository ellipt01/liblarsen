/*
 * linreg.c
 *
 *  Created on: 2014/05/19
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <linreg.h>

#include "linreg_private.h"

const int		ione  =  1;
const double	dzero =  0.;
const double	done  =  1.;
const double	dmone = -1.;

/* print an error message and exit */
void
linreg_error (const char * function_name, const char *error_msg, const char *file, const int line)
{
	fprintf (stderr, "ERROR: %s: %s:%d: %s\n", function_name, file, line, error_msg);
	exit (1);
}

linreg *
linreg_alloc (const double lambda2, const size_t n, const size_t p, const double *y, const double *x)
{
	size_t		np;
	linreg		*lreg;

	if (lambda2 < 0.) linreg_error ("lisys_alloc", "lambda2 must be >= 0.", __FILE__, __LINE__);
	if (!y) linreg_error ("lisys_alloc", "vector *y is empty.", __FILE__, __LINE__);
	if (!x) linreg_error ("lisys_alloc", "matrix *x is empty.", __FILE__, __LINE__);
	if (n <= 0) linreg_error ("lisys_alloc", "n must be > 0.", __FILE__, __LINE__);
	if (p <= 0) linreg_error ("lisys_alloc", "p must be > 0.", __FILE__, __LINE__);

	lreg = (linreg *) malloc (sizeof (linreg));

	lreg->n = n;
	lreg->p = p;
	np = lreg->n * lreg->p;

	lreg->y = (double *) malloc (lreg->n * sizeof (double));
	lreg->x = (double *) malloc (np * sizeof (double));

	dcopy_ (LINREG_CINTP (n), y, &ione, lreg->y, &ione);
	dcopy_ (LINREG_CINTP (np), x, &ione, lreg->x, &ione);

	/* By default, data is assumed to be not normalized and not standardized */
	lreg->meany = NULL;
	lreg->meanx = NULL;
	lreg->normx = NULL;
	lreg->ycentered = false;
	lreg->xcentered = false;
	lreg->xnormalized = false;

	lreg->lambda2 = lambda2;

	/* regression type */
	if (lambda2 > dlamch_ ("e")) {	// default: ridge regression
		lreg->regtype = REGULARIZATION_RIDGE;
		lreg->scale2 = 1. / (1. + lambda2);
		lreg->scale = sqrt (lreg->scale2);
	} else {	// lasso
		lreg->regtype = REGULARIZATION_LASSO;
		lreg->scale2 = 1. ;
		lreg->scale = 1.;
	}
	lreg->pen = NULL;

	return lreg;
}

void
linreg_free (linreg *lreg)
{
	if (lreg) {
		if (lreg->y) free (lreg->y);
		if (lreg->x) free (lreg->x);

		if (lreg->meany) free (lreg->meany);
		if (lreg->meanx) free (lreg->meanx);
		if (lreg->normx) free (lreg->normx);

		free (lreg);
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
		double	*xj = x + LINREG_INDEX_OF_MATRIX (0, j, size1);
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
		double	*xj = x + LINREG_INDEX_OF_MATRIX (0, j, size1);
		nrm[j] = dnrm2_ (LINREG_CINTP (size1), xj, &ione);
		alpha = 1. / nrm[j];
		dscal_ (LINREG_CINTP (size1), &alpha, xj, &ione);
	}
	return nrm;
}

/* centering lreg->y,
 * and set lreg->meany[0] = mean( y ) */
void
linreg_centering_y (linreg *lreg)
{
	lreg->meany = centering (lreg->n, 1, lreg->y);
	lreg->ycentered = true;
	return;
}

/* centering each col of lreg->x,
 * and set lreg->meanx[j] = mean( X(:,j) ) */
void
linreg_centering_x (linreg *lreg)
{
	lreg->meanx = centering (lreg->n, lreg->p, lreg->x);
	lreg->xcentered = true;
	return;
}

/* normalizing each col of lreg->x,
 * and set lreg->normx[j] = mean( X(:,j) ) */
void
linreg_normalizing_x (linreg *lreg)
{
	lreg->normx = normalizing (lreg->n, lreg->p, lreg->x);
	lreg->xnormalized = true;
	return;
}

/* standardizing lreg->x */
void
linreg_standardizing_x (linreg *lreg)
{
	linreg_centering_x (lreg);
	linreg_normalizing_x (lreg);
	return;
}

penalty *
penalty_alloc (const size_t pj, const size_t p, const double *r)
{
	penalty	*pen;
	size_t		t;

	if (!r) linreg_error ("penalty_alloc", "matrix *r is empty.", __FILE__, __LINE__);
	pen = (penalty *) malloc (sizeof (penalty));
	t = pj * p;
	pen->pj = pj;
	pen->p = p;
	pen->r = (double *) malloc (t * sizeof (double));
	dcopy_ (LINREG_CINTP (t), r, &ione, pen->r, &ione);
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

/* set lreg->pen = pen, and set lreg->scale2 = 1 / (a + b * lambda2) */
void
linreg_set_penalty (linreg *lreg, const double a, const double b, const penalty *pen)
{
	if (!pen) linreg_error ("linreg_set_penalty", "penalty *pen is empty.", __FILE__, __LINE__);
	if (lreg->p != pen->p)
		linreg_error ("linreg_set_penalty", "penalty *pen->p must be same as linreg *lreg->p.", __FILE__, __LINE__);

	if (lreg->regtype != REGULARIZATION_LASSO) {
		lreg->regtype = REGULARIZATION_USER_DEFINED;
		lreg->scale2 = 1. / (a + b * lreg->lambda2);
		lreg->scale = sqrt (lreg->scale2);
	}
	lreg->pen = pen;
}

/* is lasso? : lambda2 <= eps? */
bool
linreg_is_regtype_lasso (const linreg *lreg)
{
	return (lreg->regtype == REGULARIZATION_LASSO);
}

/* is Ridge? :
 * lambda2 > eps && default lreg->pen == NULL, i.e. default regression type? */
bool
linreg_is_regtype_ridge (const linreg *lreg)
{
	return (lreg->regtype == REGULARIZATION_RIDGE);
}

double
linreg_get_lambda2 (const linreg *lreg)
{
	return lreg->lambda2;
}

double
linreg_get_scale (const linreg *lreg)
{
	return lreg->scale;
}

double
linreg_get_scale2 (const linreg *lreg)
{
	return lreg->scale2;
}

const size_t
linreg_get_n (const linreg *lreg)
{
	return lreg->n;
}

const size_t
linreg_get_p (const linreg *lreg)
{
	return lreg->p;
}

/* return number of rows of penalty matrix: lreg->pen->pj */
const size_t
linreg_get_pj (const linreg *lreg)
{
	if (lreg->pen) return lreg->pen->pj;
	return lreg->p;
}

/* return lreg->x */
const double
*linreg_get_x (const linreg *lreg)
{
	return lreg->x;
}

/* return lreg->y */
const double
*linreg_get_y (const linreg *lreg)
{
	return lreg->y;
}

/* return lreg->pen->r */
const double
*linreg_get_penalty (const linreg *lreg)
{
	if (lreg->pen) return lreg->pen->r;
	return NULL;
}
