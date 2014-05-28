/*
 * utils.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <larsen.h>

#include "larsen_private.h"

/* print an error message and exit */
void
larsen_error (const char * function_name, const char *error_msg)
{
	fprintf (stderr, "ERROR: %s: %s\n", function_name, error_msg);
	exit (1);
}

static double
larsen_double_eps (void)
{
	char	cmach = 'e';
	return dlamch_ (&cmach);
}

static void
larsen_array_set_all (const size_t size, double *x, double val)
{
	int		i;
	for (i = 0; i < size; i++) x[i] = val;
}

larsen *
larsen_alloc (const size_t n, const size_t p, const double *y, const double *x, const double lambda1, const double lambda2)
{
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

	l->is_lasso = (lambda2 > larsen_double_eps ()) ? false : true;
	l->scale2 = (l->is_lasso) ? 1. : 1. / (1. + lambda2);
	l->scale = (l->is_lasso) ? 1. : sqrt (l->scale2);

	/* correlation */
	l->sup_c = 0.;
	l->c = NULL;

	/* active set */
	l->oper.action = ACTIVESET_ACTION_ADD;
	l->oper.index_of_A = -1;
	l->oper.column_of_X = -1;

	l->sizeA = 0;
	l->A = (int *) malloc (l->p * sizeof (int));

	/* active set */
	l->absA = 0.;
	l->u = NULL;
	l->w = NULL;

	/* solution */
	l->beta = (double *) malloc (l->p * sizeof (double));
	larsen_array_set_all (l->p, l->beta, 0.);
	l->mu = (double *) malloc (l->n * sizeof (double));
	larsen_array_set_all (l->n, l->mu, 0.);

	/* backup of solution */
	l->beta_prev = (double *) malloc (l->p * sizeof (double));
	l->mu_prev = (double *) malloc (l->n * sizeof (double));

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
		if (l->chol) free (l->chol);
		free (l);
	}
	return;
}

/* return copy of navie solution: beta_navie = l->beta */
static double *
larsen_copy_beta_navie (const larsen *l)
{
	size_t	p = l->p;
	double	*beta = (double *) malloc (p * sizeof (double));

	if (larsen_does_need_interpolation (l)) {
		// interpolation of beta
		double	stepsize = l->absA * (larsen_get_lambda1 (l, true) - l->nrm1_prev);
		dcopy_ (CINTP (p), l->beta_prev, &ione, beta, &ione);
		larsen_awpy (l, stepsize, l->w, beta);	// beta(A) += stepsize * w(A)
	} else {
		dcopy_ (CINTP (p), l->beta, &ione, beta, &ione);
	}
	return beta;
}

/* return copy of elastic net solution: beta_elnet = beta_navie / scale */
double *
larsen_copy_beta (const larsen *l, bool scale)
{
	double	*beta = larsen_copy_beta_navie (l);
	if (scale && !l->is_lasso) {
		double	alpha = 1. / l->scale;
		dscal_ (CINTP (l->p), &alpha, beta, &ione);
	}
	return beta;
}

/* return copy of navie solution: mu_navie = l->mu */
static double *
larsen_copy_mu_navie (const larsen *l)
{
	size_t	n = l->n;
	double	*mu = (double *) malloc (l->n * sizeof (double));

	if (larsen_does_need_interpolation (l)) {
		// interpolation of beta
		double	stepsize = l->absA * (larsen_get_lambda1 (l, true) - l->nrm1_prev);
		dcopy_ (CINTP (n), l->mu_prev, &ione, mu, &ione);
		daxpy_ (CINTP (n), &stepsize, l->u, &ione, mu, &ione);	// mu += stepsize * u
	} else {
		dcopy_ (CINTP (n), l->mu, &ione, mu, &ione);
	}
	return mu;
}

/* return copy of elastic net solution: mu_elnet = mu_navie / scale^2 */
double *
larsen_copy_mu (const larsen *l, bool scale)
{
	double	*mu = larsen_copy_mu_navie (l);
	if (scale && !l->is_lasso) {
		double	alpha = 1. / l->scale2;
		dscal_ (CINTP (l->n), &alpha, mu, &ione);
	}
	return mu;
}

/* increment l->lambda1 */
void
larsen_set_lambda1 (larsen *l, double t)
{
	l->lambda1 = t;
	return;
}

/* return l->lambda1 */
double
larsen_get_lambda1 (const larsen *l, bool scaling)
{
	double	lambda1 = l->lambda1;
	if (scaling && !l->is_lasso) lambda1 *= l->scale;
	return lambda1;
}
