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
larsen_alloc (const linsys *lsys, const double lambda1)
{
	larsen	*l;

	if (!lsys) linsys_error ("larsen_alloc", "linsys *lsys is empty.", __FILE__, __LINE__);
	if (lambda1 < 0) linsys_error ("larsen_alloc", "lambda1 must be >= 0.", __FILE__, __LINE__);

	if (!lsys->ycentered)
		linsys_error ("larsen_alloc", "for LARSE-EN, vector *y must be centered.\nplease use linsys_centering_y()", __FILE__, __LINE__);
	if (!lsys->xcentered || !lsys->xnormalized)
		linsys_error ("larsen_alloc", "for LARSE-EN, matrix *x must be standardized.\nplease use linsys_{centering,normalizing}_x()", __FILE__, __LINE__);

	l = (larsen *) malloc (sizeof (larsen));

	l->stop_loop = false;

	/* register linear system of regression equations */
	l->lsys = lsys;

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
	l->A = (int *) malloc (l->lsys->p * sizeof (int));

	/* active set */
	l->absA = 0.;
	l->u = NULL;
	l->w = NULL;

	/* solution */
	l->beta = (double *) malloc (l->lsys->p * sizeof (double));
	array_set_all (l->lsys->p, l->beta, 0.);
	l->mu = (double *) malloc (l->lsys->n * sizeof (double));
	array_set_all (l->lsys->n, l->mu, 0.);

	/* backup of solution */
	l->beta_prev = (double *) malloc (l->lsys->p * sizeof (double));
	l->mu_prev = (double *) malloc (l->lsys->n * sizeof (double));

	/* interpolation */
	l->is_interped = false;
	l->stepsize_intr = 0.;
	l->beta_intr = (double *) malloc (l->lsys->p * sizeof (double));
	l->mu_intr = (double *) malloc (l->lsys->n * sizeof (double));

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
	size_t	p = linsys_get_p (l->lsys);
	double	*beta = (double *) malloc (p * sizeof (double));
	if (!l->is_interped) dcopy_ (LINSYS_CINTP (p), l->beta, &ione, beta, &ione);
	else dcopy_ (LINSYS_CINTP (p), l->beta_intr, &ione, beta, &ione);
	return beta;
}

/* return copy of beta
 * if scaling == true && !lasso, return copy of l->beta / scale */
double *
larsen_copy_beta (const larsen *l, bool scaling)
{
	size_t	p = linsys_get_p (l->lsys);
	double	scale = linsys_get_scale (l->lsys);
	double	*beta = larsen_copy_beta_navie (l);
	if (scaling && !linsys_is_regtype_lasso (l->lsys)) {
		double	alpha = 1. / scale;
		dscal_ (LINSYS_CINTP (p), &alpha, beta, &ione);
	}
	return beta;
}

/* return copy of l->mu */
static double *
larsen_copy_mu_navie (const larsen *l)
{
	size_t	n = linsys_get_n (l->lsys);
	double	*mu = (double *) malloc (n * sizeof (double));
	if (!l->is_interped) dcopy_ (LINSYS_CINTP (n), l->mu, &ione, mu, &ione);
	else dcopy_ (LINSYS_CINTP (n), l->mu_intr, &ione, mu, &ione);
	return mu;
}

/* return copy of mu
 * if scaling == true && !lasso, return copy of l->mu / scale^2 */
double *
larsen_copy_mu (const larsen *l, bool scaling)
{
	size_t	n = linsys_get_n (l->lsys);
	double	scale2 = linsys_get_scale2 (l->lsys);
	double	*mu = larsen_copy_mu_navie (l);
	if (scaling && !linsys_is_regtype_lasso (l->lsys)) {
		double	alpha = 1. / scale2;
		dscal_ (LINSYS_CINTP (n), &alpha, mu, &ione);
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
	if (scaling && !linsys_is_regtype_lasso (l->lsys)) lambda1 *= linsys_get_scale (l->lsys);
	return lambda1;
}
