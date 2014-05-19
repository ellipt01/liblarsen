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

#include "linsys_private.h"

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
larsen_alloc (const linsys *sys, const double lambda1, const double lambda2)
{
	larsen	*l;

	if (!sys) linsys_error ("larsen_alloc", "linsys *sys is empty.", __FILE__, __LINE__);
	if (lambda1 < 0 || lambda2 < 0) {
		char	msg[80];
		sprintf (msg, "ERROR: lambda1(= %f), lambda2(= %f) must be >= 0.", lambda1, lambda2);
		linsys_error ("larsen_alloc", msg, __FILE__, __LINE__);
	}

	l = (larsen *) malloc (sizeof (larsen));

	l->sys = sys;
	if (!l->sys->y_centerdized)
		linsys_error ("larsen_alloc", "vector *y must be centerdized.", __FILE__, __LINE__);
	if (!l->sys->x_centerdized || !l->sys->x_centerdized)
		linsys_error ("larsen_alloc", "matrix *x must be normalized.", __FILE__, __LINE__);

	l->stop_loop = false;

	l->lambda1 = lambda1;
	l->lambda2 = lambda2;

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
	l->A = (int *) malloc (l->sys->p * sizeof (int));

	/* active set */
	l->absA = 0.;
	l->u = NULL;
	l->w = NULL;

	/* solution */
	l->beta = (double *) malloc (l->sys->p * sizeof (double));
	larsen_array_set_all (l->sys->p, l->beta, 0.);
	l->mu = (double *) malloc (l->sys->n * sizeof (double));
	larsen_array_set_all (l->sys->n, l->mu, 0.);

	/* backup of solution */
	l->beta_prev = (double *) malloc (l->sys->p * sizeof (double));
	l->mu_prev = (double *) malloc (l->sys->n * sizeof (double));

	/* interpolation */
	l->interp = false;
	l->stepsize_intr = 0.;
	l->beta_intr = (double *) malloc (l->sys->p * sizeof (double));
	l->mu_intr = (double *) malloc (l->sys->n * sizeof (double));

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
larsen_copy_beta_navie (const larsen *l)
{
	double	*beta = (double *) malloc (l->sys->p * sizeof (double));
	if (!l->interp) dcopy_ (LINSYS_CINTP (l->sys->p), l->beta, &ione, beta, &ione);
	else dcopy_ (LINSYS_CINTP (l->sys->p), l->beta_intr, &ione, beta, &ione);
	return beta;
}

/* return copy of elastic net solution: beta_elnet = beta_navie / scale */
double *
larsen_copy_beta_elasticnet (const larsen *l)
{
	double	*beta = larsen_copy_beta_navie (l);
	if (!l->is_lasso) {
		double	alpha = 1. / l->scale;
		dscal_ (LINSYS_CINTP (l->sys->p), &alpha, beta, &ione);
	}
	return beta;
}

/* return copy of navie solution: mu_navie = l->mu */
double *
larsen_copy_mu_navie (const larsen *l)
{
	double	*mu = (double *) malloc (l->sys->n * sizeof (double));
	if (!l->interp) dcopy_ (LINSYS_CINTP (l->sys->n), l->mu, &ione, mu, &ione);
	else dcopy_ (LINSYS_CINTP (l->sys->n), l->mu_intr, &ione, mu, &ione);
	return mu;
}

/* return copy of elastic net solution: mu_elnet = mu_navie / scale^2 */
double *
larsen_copy_mu_elasticnet (const larsen *l)
{
	double	*mu = larsen_copy_mu_navie (l);
	if (!l->is_lasso) {
		double	alpha = 1. / l->scale2;
		dscal_ (LINSYS_CINTP (l->sys->n), &alpha, mu, &ione);
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
