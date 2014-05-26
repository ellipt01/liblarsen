/*
 * larsen.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdlib.h>
#include <math.h>
#include <larsen.h>

#include "larsen_private.h"

/* Calculate correlation between residual (y - X * beta) and predictors X' */
static void
update_correlations (larsen *l)
{
	size_t			n = l->lreg->n;
	size_t			p = l->lreg->p;
	double			scale = l->lreg->scale;
	const double	*x = l->lreg->x;
	const double	*y = l->lreg->y;

	/*
	 *  c = Z' * (b - Z * beta),  b = [y; 0], Z = scale * [X; sqrt(lambda2) * J]
	 *  if lambda2 > 0 (scale != 1),
	 *  c = scale * X' * (y - mu) - scale^2 * lambda2 * J' * J * beta
	 */

	// c = scale * X' * (y - mu)
	double	*r = (double *) malloc (n * sizeof (double));
	dcopy_ (LINREG_CINTP (n), y, &ione, r, &ione);
	daxpy_ (LINREG_CINTP (n), &dmone, l->mu, &ione, r, &ione);	// r = - mu + y

	if (l->c) free (l->c);
	l->c = (double *) malloc (p * sizeof (double));
	dgemv_ ("T", LINREG_CINTP (n), LINREG_CINTP (p), &scale, x, LINREG_CINTP (n), r, &ione, &dzero, l->c, &ione);
	free (r);

	/*** c -= scale2 * lambda2 * J' * J * beta ***/
	if (!larsen_is_regtype_lasso (l)) {	// lambda2 > 0
		double		alpha = - l->lreg->lambda2 * l->lreg->scale2;

		if (larsen_is_regtype_ridge (l)) {	// Ridge
			/*** c -= scale2 * lambda2 * E' * E * beta ***/
			daxpy_ (LINREG_CINTP (p), &alpha, l->beta, &ione, l->c, &ione);
		} else {
			/*** c -= scale2 * lambda2 * J' * J * beta ***/
			size_t			pj = l->lreg->pen->pj;
			const double	*jr = l->lreg->pen->r;
			double			*jb = (double *) malloc (pj * sizeof (double));
			double			*jtjb = (double *) malloc (p * sizeof (double));

			// J * beta
			dgemv_ ("N", LINREG_CINTP (pj), LINREG_CINTP (p), &done, jr, LINREG_CINTP (pj), l->beta, &ione, &dzero, jb, &ione);
			// J' * (J * beta)
			dgemv_ ("T", LINREG_CINTP (pj), LINREG_CINTP (p), &done, jr, LINREG_CINTP (pj), jb, &ione, &dzero, jtjb, &ione);
			free (jb);
			// c -= scale2 * lambda2 * J' * J * beta
			daxpy_ (LINREG_CINTP (p), &alpha, jtjb, &ione, l->c, &ione);
			free (jtjb);
		}
	}

	{
		int		maxidx = idamax_ (LINREG_CINTP (p), l->c, &ione) - 1;	// differ of fortran and C
		if (l->sizeA == 0) {
			l->oper.index_of_A = 0;
			l->oper.column_of_X = maxidx;
		}
		l->sup_c = fabs (l->c[maxidx]);
	}

	return;
}

/* update beta: beta += stepsize * l->w */
void
update_beta (const larsen *l, const double stepsize, double *beta)
{
	larsen_axapy (l, stepsize, l->w, beta);	// beta(A) += stepsize * w(A)
	return;
}

/* update mu: mu += stepsize * l->u */
void
update_mu (const larsen *l, const double stepsize, double *mu)
{
	int		n = (int) l->lreg->n;
	daxpy_ (&n, &stepsize, l->u, &ione, mu, &ione);	// mu += stepsize * u
	return;
}

/* Update beta and mu. beta += stepsize * w, mu += stepsze * u */
static void
update_solutions (larsen *l)
{
	size_t		n = l->lreg->n;
	size_t		p = l->lreg->p;

	// backup previous results
	l->nrm1_prev = l->nrm1;
	dcopy_ (LINREG_CINTP (p), l->beta, &ione, l->beta_prev, &ione);
	dcopy_ (LINREG_CINTP (n), l->mu, &ione, l->mu_prev, &ione);

	// update beta and mu
	update_beta (l, l->stepsize, l->beta);
	update_mu (l, l->stepsize, l->mu);
	l->nrm1 = dasum_ (LINREG_CINTP (p), l->beta, &ione);

	return;
}

static void
update_stop_loop_flag (larsen *l)
{
	size_t	n = l->lreg->n;
	size_t	p = l->lreg->p;
	int		size = l->sizeA;
	int		m = (larsen_is_regtype_lasso (l)) ? (int) LINREG_MIN (n - 1, p) : (int) p;
	if (l->oper.action == ACTIVESET_ACTION_DROP) size--;
	l->stop_loop = (size >= m) ? true : false;
	if (!l->stop_loop) l->stop_loop = (l->oper.column_of_X == -1);
	return;
}

/* Progress one step of the LARS-EN algorithm
 * add / remove one variable (assigned by l->oper.column)
 * to / from the active set and update equiangular vector, stepsize,
 * and solution vectors. Also, stop_loop flag is modified.
 *
 * 1. Update correlation between residual (y - mu) and variables (X').
 *    In the case of A = {} (on first iteration or restarted), a variable
 *    which has largest correlation is selected as the active set item
 *    (l->oper.column_of_X is set to its index in X).
 *
 * 2. Add / remove one variable assigned by l->oper.column_of_X to / from the active set.
 *    When l->oper.action == ACTIVESET_ACTION_ADD (add a new item to active set),
 *    the value of l->oper.column_of_X is stored in l->A[l->oper.index_of_A]
 *    and when l->oper.action == ACTIVESET_ACTION_DROP (remove a item from active set)
 *    l->A[l->oper.index_of_A] ( == l->oper.column_of_X) is removed.
 *
 * 3. Update equi-angular vector and its relevant.
 *
 * 4. Update step size and l->oper (which specify the next operation to the active set).
 *    l->oper.action is updated according to whether gamma_hat < gamma_tilde or not.
 *    Case gamma_tilde < gamma_hat:
 *      l->oper.index_of_A  = argminplus_i ( gamma_tilde( A(i) ) )
 *      l->oper.column_of_X = argminplus_j ( gamma_tilde( X(j) ) )
 *    Case gamma_tilde > gamma_hat:
 *      l->oper.index_of_A  = argminplus_i ( gamma_hat( Ac(i) ) )
 *      l->oper.column_of_X = argminplus_j ( gamma_hat( X (j) ) )
 *
 * 5. Update solutions (beta and mu).
 *
 * 6. Update stop_loop flag. if it sets to true, loop of the regression should be terminated.
 */
bool
larsen_regression_step (larsen *l)
{
	l->stop_loop = true;

	update_correlations (l);

	if (!update_activeset (l)) return false;

	if (!update_equiangular (l)) return false;

	if (!update_stepsize (l)) return false;

	update_solutions (l);

	update_stop_loop_flag (l);

	return true;
}

/*
 * check whether interpolation does need.
 * if nrm1_prev <= scale * lambda < nrm1, need interpolation
 */
bool
larsen_does_need_interpolation (const larsen *l)
{
	double	lambda1 = larsen_get_lambda1 (l, true);
	return (l->nrm1_prev <= lambda1 && lambda1 < l->nrm1);
}
