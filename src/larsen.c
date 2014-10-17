/*
 * utils.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <larsen.h>

#include "private.h"

extern void	larsen_awpy (const larsen *l, double alpha, double *w, double *y);

/* print an error message and exit */
void
larsen_error (const char * function_name, const char *error_msg)
{
	fprintf (stderr, "ERROR: %s: %s\n", function_name, error_msg);
	exit (1);
}

static larsen *
larsen_alloc (void)
{
	larsen	*l = (larsen *) malloc (sizeof (larsen));
	l->c = NULL;
	l->beta = NULL;
	l->mu = NULL;
	l->beta_prev = NULL;
	l->mu_prev = NULL;
	return l;
}

larsen *
larsen_new (const linregmodel *lreg, const double *x, const double lambda1)
{
	larsen	*l;

	if (!lreg) error_and_exit ("larsen_new", "linregmodel is empty.", __FILE__, __LINE__);
	if (mm_real_is_sparse (lreg->x)) error_and_exit ("larsen_new", "sparse matrix is not supported.", __FILE__, __LINE__);
	if (mm_real_is_symmetric (lreg->x)) error_and_exit ("larsen_new", "symmetric matrix is not supported.", __FILE__, __LINE__);
	if (lambda1 < 0) error_and_exit ("larsen_new", "lambda1 must be >= 0.", __FILE__, __LINE__);

	l = larsen_alloc ();

	l->lreg = lreg;

	l->stop_loop = false;

	l->lambda1 = lambda1;

	l->scale2 = (l->lreg->is_regtype_lasso) ? 1. : 1. / (1. + lreg->lambda2);
	l->scale = (l->lreg->is_regtype_lasso) ? 1. : sqrt (l->scale2);

	/* correlation */
	l->sup_c = 0.;
	l->c = NULL;

	/* active set */
	l->oper.action = ACTIVESET_ACTION_ADD;
	l->oper.index_of_A = -1;
	l->oper.column_of_X = -1;

	l->sizeA = 0;
	l->A = (int *) malloc (l->lreg->x->n * sizeof (int));

	/* active set */
	l->absA = 0.;
	l->u = NULL;
	l->w = NULL;

	/* solution */
	l->beta = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, l->lreg->x->n, 1, l->lreg->x->n);
	l->beta->data = (double *) malloc (l->beta->nz * sizeof (double));
	mm_real_set_all (l->beta, 0.);
	l->mu = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, l->lreg->y->m, 1, l->lreg->y->nz);
	l->mu->data = (double *) malloc (l->mu->nz * sizeof (double));
	mm_real_set_all (l->mu, 0.);

	/* backup of solution */
	l->beta_prev = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, l->lreg->x->n, 1, l->lreg->x->n);
	l->beta_prev->data = (double *) malloc (l->beta_prev->nz * sizeof (double));
	l->mu_prev = mm_real_new (MM_REAL_DENSE, MM_REAL_GENERAL, l->lreg->y->m, 1, l->lreg->y->nz);
	l->mu_prev->data = (double *) malloc (l->mu_prev->nz * sizeof (double));

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
		if (l->beta) mm_real_free (l->beta);
		if (l->mu) free (l->mu);
		if (l->beta_prev) mm_real_free (l->beta_prev);
		if (l->mu_prev) free (l->mu_prev);
		if (l->chol) free (l->chol);
		free (l);
	}
	return;
}

/* If l->nrm1_prev < lambda1 < l->nrm1, need interpolation */
static bool
larsen_does_need_interpolation (const larsen *l)
{
	double	lambda1 = larsen_get_lambda1 (l, true);
	return (l->nrm1_prev <= lambda1 && lambda1 < l->nrm1);
}

/* return copy of navie solution: beta_navie = l->beta */
static double *
larsen_copy_beta_navie (const larsen *l)
{
	double	*beta = (double *) malloc (l->beta->nz * sizeof (double));

	if (larsen_does_need_interpolation (l)) {
		// interpolation of beta
		double	stepsize = l->absA * (larsen_get_lambda1 (l, true) - l->nrm1_prev);
		dcopy_ (&l->beta_prev->nz, l->beta_prev->data, &ione, beta, &ione);
		larsen_awpy (l, stepsize, l->w, beta);	// beta(A) += stepsize * w(A)
	} else {
		dcopy_ (&l->beta->nz, l->beta->data, &ione, beta, &ione);
	}
	return beta;
}

/* return copy of elastic net solution: beta_elnet = beta_navie / scale */
double *
larsen_copy_beta (const larsen *l, bool scale)
{
	double	*beta = larsen_copy_beta_navie (l);
	if (scale && !l->lreg->is_regtype_lasso) {
		double	alpha = 1. / l->scale;
		dscal_ (&l->beta->nz, &alpha, beta, &ione);
	}
	return beta;
}

/* return copy of navie solution: mu_navie = l->mu */
static double *
larsen_copy_mu_navie (const larsen *l)
{
//	size_t	n = l->n;
	double	*mu = (double *) malloc (l->mu->nz * sizeof (double));

	if (larsen_does_need_interpolation (l)) {
		// interpolation of beta
		double	stepsize = l->absA * (larsen_get_lambda1 (l, true) - l->nrm1_prev);
		dcopy_ (&l->mu_prev->nz, l->mu_prev->data, &ione, mu, &ione);
		daxpy_ (&l->mu->nz, &stepsize, l->u, &ione, mu, &ione);	// mu += stepsize * u
	} else {
		dcopy_ (&l->mu->nz, l->mu->data, &ione, mu, &ione);
	}
	return mu;
}

/* return copy of elastic net solution: mu_elnet = mu_navie / scale^2 */
double *
larsen_copy_mu (const larsen *l, bool scale)
{
	double	*mu = larsen_copy_mu_navie (l);
	if (scale && !l->lreg->is_regtype_lasso) {
		double	alpha = 1. / l->scale2;
		dscal_ (&l->mu->nz, &alpha, mu, &ione);
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
	if (scaling && !l->lreg->is_regtype_lasso) lambda1 *= l->scale;
	return lambda1;
}
