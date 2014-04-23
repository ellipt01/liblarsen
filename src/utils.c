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

/* return copy of navie solution: beta_navie = l->beta */
double *
larsen_copy_beta_navie (larsen *l)
{
	double	*beta = (double *) malloc (l->p * sizeof (double));
	if (!l->interp) cblas_dcopy (l->p, l->beta, 1, beta, 1);
	else cblas_dcopy (l->p, l->beta_intr, 1, beta, 1);
	return beta;
}

/* return copy of elastic net solution: beta_elnet = beta_navie / scale */
double *
larsen_copy_beta_elasticnet (larsen *l)
{
	double	*beta = larsen_copy_beta_navie (l);
	if (l->is_elnet) cblas_dscal (l->p, 1. / l->scale, beta, 1);
	return beta;
}

/* return copy of navie solution: mu_navie = l->mu */
double *
larsen_copy_mu_navie (larsen *l)
{
	double	*mu = (double *) malloc (l->n * sizeof (double));
	if (!l->interp) cblas_dcopy (l->n, l->mu, 1, mu, 1);
	else cblas_dcopy (l->n, l->mu_intr, 1, mu, 1);
	return mu;
}

/* return copy of elastic net solution: mu_elnet = mu_navie / scale^2 */
double *
larsen_copy_mu_elasticnet (larsen *l)
{
	double	*mu = larsen_copy_mu_navie (l);
	if (l->is_elnet) cblas_dscal (l->n, 1. / l->scale2, mu, 1);
	return mu;
}

/* increment l->lambda1 */
void
larsen_set_lambda1 (larsen *l, double t)
{
	l->lambda1 = t;
	return;
}
