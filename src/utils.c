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
	l->beta = (double *) malloc (l->lreg->p * sizeof (double));
	array_set_all (l->lreg->p, l->beta, 0.);
	l->mu = (double *) malloc (l->lreg->n * sizeof (double));
	array_set_all (l->lreg->n, l->mu, 0.);

	/* backup of solution */
	l->beta_prev = (double *) malloc (l->lreg->p * sizeof (double));
	l->mu_prev = (double *) malloc (l->lreg->n * sizeof (double));

	/* interpolation */
	l->is_interped = false;
	l->stepsize_intr = 0.;
	l->beta_intr = (double *) malloc (l->lreg->p * sizeof (double));
	l->mu_intr = (double *) malloc (l->lreg->n * sizeof (double));

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
		if (l->beta_intr) free (l->beta_intr);
		if (l->mu_intr) free (l->mu_intr);
		if (l->chol) free (l->chol);
		free (l);
	}
	return;
}

/* return copy of l->beta */
static double *
larsen_copy_beta_navie (const larsen *l)
{
	size_t	p = linreg_get_p (l->lreg);
	double	*beta = (double *) malloc (p * sizeof (double));
	if (l->is_interped) dcopy_ (LINREG_CINTP (p), l->beta_intr, &ione, beta, &ione);
	else dcopy_ (LINREG_CINTP (p), l->beta, &ione, beta, &ione);
	return beta;
}

/* return copy of beta
 * if scaling == true && !lasso, return copy of l->beta / scale */
double *
larsen_copy_beta (const larsen *l, bool scaling)
{
	size_t	p = linreg_get_p (l->lreg);
	double	scale = linreg_get_scale (l->lreg);
	double	*beta = larsen_copy_beta_navie (l);
	if (scaling && !linreg_is_regtype_lasso (l->lreg)) {
		double	alpha = 1. / scale;
		dscal_ (LINREG_CINTP (p), &alpha, beta, &ione);
	}
	return beta;
}

/* return copy of l->mu */
static double *
larsen_copy_mu_navie (const larsen *l)
{
	size_t	n = linreg_get_n (l->lreg);
	double	*mu = (double *) malloc (n * sizeof (double));
	if (l->is_interped) dcopy_ (LINREG_CINTP (n), l->mu_intr, &ione, mu, &ione);
	else dcopy_ (LINREG_CINTP (n), l->mu, &ione, mu, &ione);
	return mu;
}

/* return copy of mu
 * if scaling == true && !lasso, return copy of l->mu / scale^2 */
double *
larsen_copy_mu (const larsen *l, bool scaling)
{
	size_t	n = linreg_get_n (l->lreg);
	double	scale2 = linreg_get_scale2 (l->lreg);
	double	*mu = larsen_copy_mu_navie (l);
	if (scaling && !linreg_is_regtype_lasso (l->lreg)) {
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
	if (scaling && !linreg_is_regtype_lasso (l->lreg)) lambda1 *= linreg_get_scale (l->lreg);
	return lambda1;
}
