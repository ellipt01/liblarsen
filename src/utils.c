/*
 * utils.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdlib.h>
#include <math.h>
#include <larsen.h>

larsen *
larsen_alloc (size_t n, size_t p, const double *y, const double *x, double lambda1, double lambda2)
{
	int		i;
	larsen	*l;

	if (!y) {
		fprintf (stderr, "ERROR: vector *y is empty.\n");
		return NULL;
	}
	if (!x) {
		fprintf (stderr, "ERROR: matrix *x is empty.\n");
		return NULL;
	}
	if (lambda1 < 0 || lambda2 < 0) {
		fprintf (stderr, "ERROR: lambda1(= %f), lambda2(= %f) must be >= 0.\n", lambda1, lambda2);
		return NULL;
	}

	l = (larsen *) malloc (sizeof (larsen));

	l->n = n;
	l->p = p;

	l->stop_loop = false;

	l->lambda1 = lambda1;
	l->lambda2 = lambda2;

	/* vector and matrix view of original y and x */
	l->x = x;
	l->y = y;

	l->is_elnet = (lambda2 > DBL_EPSILON);
	l->scale2 = (l->is_elnet) ? 1. / (1. + lambda2) : 1.;
	l->scale = (l->is_elnet) ? sqrt (l->scale2) : 1.;

	/* correlation */
	l->sup_c = 0.;
	l->c = NULL;

	/* active set */
	l->oper.action = ACTIVESET_ACTION_ADD;
	l->oper.index_of_A = -1;
	l->oper.column_of_X = -1;

	l->sizeA = 0;
	l->A = (int *) malloc (p * sizeof (int));

	/* active set */
	l->absA = 0.;
	l->u = NULL;
	l->w = NULL;

	/* solution */
	l->beta = (double *) malloc (p * sizeof (double));
	l->mu = (double *) malloc (n * sizeof (double));
	for (i = 0; i < p; i++) l->beta[i] = 0.;
	for (i = 0; i < n; i++) l->mu[i] = 0.;

	/* backup of solution */
	l->beta_prev = (double *) malloc (p * sizeof (double));
	l->mu_prev = (double *) malloc (n * sizeof (double));

	/* interpolation */
	l->interp = false;
	l->stepsize_intr = 0.;
	l->beta_intr = (double *) malloc (p * sizeof (double));
	l->mu_intr = (double *) malloc (n * sizeof (double));

	/* cholesky factorization */
	l->chol = NULL;

	return l;
}

void
larsen_free (larsen *l)
{
	if (l) {
		if (l->c) free (l->c);

		if (l->A) free (l->A);
//		if (l->Ac) free (l->Ac);

		if (l->u) free (l->u);
		if (l->w) free (l->w);

		if (l->beta) free (l->beta);
		if (l->mu) free (l->mu);

		if (l->beta_prev) free (l->beta_prev);
		if (l->mu_prev) free (l->mu_prev);

		if (l->beta_intr) free (l->beta_intr);
		if (l->mu_intr) free (l->mu_intr);

		if (l->chol) free (l->chol);

		free (l);
	}
	return;
}
