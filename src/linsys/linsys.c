/*
 * linsys.c
 *
 *  Created on: 2014/05/19
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
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
linsys_alloc (const size_t n, const size_t p, const double *y, const double *x, const double *meany, const double *meanx, const double *normx)
{
	size_t		np;
	linsys		*ls;

	if (!y) linsys_error ("lisys_alloc", "vector *y is empty.", __FILE__, __LINE__);
	if (!x) linsys_error ("lisys_alloc", "matrix *x is empty.", __FILE__, __LINE__);
	if (n <= 0) linsys_error ("lisys_alloc", "n must be > 0.", __FILE__, __LINE__);
	if (p <= 0) linsys_error ("lisys_alloc", "p must be > 0.", __FILE__, __LINE__);

	ls = (linsys *) malloc (sizeof (linsys));

	ls->n = n;
	ls->p = p;
	np = ls->n * ls->p;

	ls->y = (double *) malloc (ls->n * sizeof (double));
	ls->x = (double *) malloc (np * sizeof (double));

	dcopy_ (LINSYS_CINTP (ls->n), y, &ione, ls->y, &ione);
	dcopy_ (LINSYS_CINTP (np), x, &ione, ls->x, &ione);

	ls->meany = NULL;
	if (meany) {
		ls->meany = (double *) malloc (sizeof (double));
		ls->meany[0] = meany[0];
	}
	ls->y_centerdized = (ls->meany != NULL);


	ls->meanx = NULL;
	if (meanx) {
		ls->meanx = (double *) malloc (ls->p * sizeof (double));
		dcopy_ (LINSYS_CINTP (ls->p), meanx, &ione, ls->meanx, &ione);
	}
	ls->x_centerdized = (ls->meanx != NULL);

	ls->normx = NULL;
	if (normx) {
		ls->normx = (double *) malloc (ls->p * sizeof (double));
		dcopy_ (LINSYS_CINTP (ls->p), normx, &ione, ls->normx, &ione);
	}
	ls->x_normalized = (ls->normx != NULL);

	ls->pen = NULL;

	return ls;
}

void
linsys_free (linsys *ls)
{
	if (ls) {
		if (ls->y) free (ls->y);
		if (ls->x) free (ls->x);

		if (ls->meany) free (ls->meany);
		if (ls->meanx) free (ls->meanx);
		if (ls->normx) free (ls->normx);

		free (ls);
	}
	return;
}

/* Centering each column of matrix:
 * x(:, j) -> x(:, j) - mean(x(:, j)) */
double *
linsys_centering (const size_t size1, const size_t size2, double *x)
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
double *
linsys_normalizing (const size_t size1, const size_t size2, double *x)
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

penalty *
penalty_alloc (const size_t p1, const size_t p, const double a, const double b, const double *r)
{
	penalty	*pen = (penalty *) malloc (sizeof (penalty));
	size_t		t = p1 * p;
	pen->p1 = p1;
	pen->a = a;
	pen->b = b;
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

void
linsys_set_penalty (linsys *ls, const penalty *pen)
{
	ls->pen = pen;
	return;
}
