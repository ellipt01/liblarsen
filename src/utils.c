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

static void
array_set_all (const size_t size, double *x, double val)
{
	int		i;
	for (i = 0; i < size; i++) x[i] = val;
}

larsen *
larsen_alloc (const linreg *lreg, const double lambda1)
{
	larsen	*l;

	if (!lreg) linreg_error ("larsen_alloc", "linreg *lreg is empty.", __FILE__, __LINE__);
	if (lambda1 < 0) linreg_error ("larsen_alloc", "lambda1 must be >= 0.", __FILE__, __LINE__);

	if (!lreg->ycentered)
		linreg_error ("larsen_alloc", "for LARSE-EN, vector *y must be centered.\nplease use linreg_centering_y()", __FILE__, __LINE__);
	if (!lreg->xcentered || !lreg->xnormalized)
		linreg_error ("larsen_alloc", "for LARSE-EN, matrix *x must be standardized.\nplease use linreg_{centering,normalizing}_x()", __FILE__, __LINE__);

	l = (larsen *) malloc (sizeof (larsen));

	l->stop_loop = false;

	/* register linear system of regression equations */
	l->lreg = lreg;

	/*L1 threshold */
	l->lambda1 = lambda1;

	/* correlation */
	l->sup_c = 0.;
	l->c = NULL;

	/* active set */
	l->oper.action = ACTIVESET_ACTION_ADD;
	l->oper.index_of_A = -1;
	l->oper.column_of_X = -1;

	l->sizeA = 0;
	l->A = (int *) malloc (l->lreg->p * sizeof (int));

	/* active set */
	l->absA = 0.;
	l->u = NULL;
	l->w = NULL;

	/* solution */
	l->nrm1 = 0.;
	l->beta = (double *) malloc (l->lreg->p * sizeof (double));
	array_set_all (l->lreg->p, l->beta, 0.);
	l->mu = (double *) malloc (l->lreg->n * sizeof (double));
	array_set_all (l->lreg->n, l->mu, 0.);

	/* backup of solution */
	l->nrm1_prev = 0.;
	l->beta_prev = (double *) malloc (l->lreg->p * sizeof (double));
	l->mu_prev = (double *) malloc (l->lreg->n * sizeof (double));

	/* Cholesky factorization of Z(:,A)' * Z(:,A) */
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

/* return copy of l->beta */
static double *
larsen_copy_beta_navie (const larsen *l)
{
	size_t	p = l->lreg->p;
	double	*beta = (double *) malloc (p * sizeof (double));

	if (larsen_does_need_interpolation (l)) {
		// interpolation of beta
		double	stepsize = l->absA * (larsen_get_lambda1 (l, true) - l->nrm1_prev);
		dcopy_ (LINREG_CINTP (p), l->beta_prev, &ione, beta, &ione);
		update_beta (l, stepsize, beta);
	} else {
		dcopy_ (LINREG_CINTP (p), l->beta, &ione, beta, &ione);
	}

	return beta;
}

/* return copy of beta
 * if scaling == true && !lasso, return copy of l->beta / scale */
double *
larsen_copy_beta (const larsen *l, bool scaling)
{
	size_t	p = l->lreg->p;
	double	scale = l->lreg->scale;
	double	*beta = larsen_copy_beta_navie (l);
	if (scaling && !larsen_is_regtype_lasso (l)) {
		double	alpha = 1. / scale;
		dscal_ (LINREG_CINTP (p), &alpha, beta, &ione);
	}
	return beta;
}

/* return copy of l->mu */
static double *
larsen_copy_mu_navie (const larsen *l)
{
	size_t	n = l->lreg->n;
	double	*mu = (double *) malloc (n * sizeof (double));

	if (larsen_does_need_interpolation (l)) {
		// interpolation of mu
		double	stepsize = l->absA * (larsen_get_lambda1 (l, true) - l->nrm1_prev);
		dcopy_ (LINREG_CINTP (n), l->mu_prev, &ione, mu, &ione);
		update_mu (l, stepsize, mu);
	} else {
		dcopy_ (LINREG_CINTP (n), l->mu, &ione, mu, &ione);
	}

	return mu;
}

/* return copy of mu
 * if scaling == true && !lasso, return copy of l->mu / scale^2 */
double *
larsen_copy_mu (const larsen *l, bool scaling)
{
	size_t	n = l->lreg->n;
	double	scale2 = l->lreg->scale2;
	double	*mu = larsen_copy_mu_navie (l);
	if (scaling && !larsen_is_regtype_lasso (l)) {
		double	alpha = 1. / scale2;
		dscal_ (LINREG_CINTP (n), &alpha, mu, &ione);
	}
	return mu;
}

/* set l->lambda1 = t */
void
larsen_set_lambda1 (larsen *l, double t)
{
	l->lambda1 = t;
	return;
}

/* return l->lambda1
 * if scaling == true && !lasso, return scale * l->lambda1 */
double
larsen_get_lambda1 (const larsen *l, bool scaling)
{
	double	lambda1 = l->lambda1;
	if (scaling && !larsen_is_regtype_lasso (l)) lambda1 *= l->lreg->scale;
	return lambda1;
}

/* lambda2 <= eps, regression type = lasso */
bool
larsen_is_regtype_lasso (const larsen *l)
{
	return (l->lreg->pentype == NO_PENALTY);
}

/* lambda2 > eps && l->lreg->pen == NULL, regression type = Ridge */
bool
larsen_is_regtype_ridge (const larsen *l)
{
	if (larsen_is_regtype_lasso (l)) return false;
	return (l->lreg->pentype == PENALTY_RIDGE);
}
