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
	size_t			n = linsys_get_n (l->lsys);
	size_t			p = linsys_get_p (l->lsys);
	double			scale = linsys_get_scale (l->lsys);
	const double	*x = linsys_get_x (l->lsys);
	const double	*y = linsys_get_y (l->lsys);

	/*
	 *  c = Z' * (b - Z * beta),  b = [y; 0], Z = scale * [X; sqrt(lambda2) * J]
	 *  if lambda2 > 0 (scale != 1),
	 *  c = scale * X' * (y - mu) - scale^2 * lambda2 * J' * J * beta
	 */

	// c = scale * X' * (y - mu)
	double	*r = (double *) malloc (n * sizeof (double));
	dcopy_ (LINSYS_CINTP (n), y, &ione, r, &ione);
	daxpy_ (LINSYS_CINTP (n), &dmone, l->mu, &ione, r, &ione);	// r = - mu + y

	if (l->c) free (l->c);
	l->c = (double *) malloc (p * sizeof (double));
	dgemv_ ("T", LINSYS_CINTP (n), LINSYS_CINTP (p), &scale, x, LINSYS_CINTP (n), r, &ione, &dzero, l->c, &ione);
	free (r);

	/*** c -= scale2 * lambda2 * J' * J * beta ***/
	if (!linsys_is_regtype_lasso (l->lsys)) {	// lambda2 > 0
		double		lambda2 = linsys_get_lambda2 (l->lsys);
		double		scale2 = linsys_get_scale2 (l->lsys);
		double		alpha = - lambda2 * scale2;

		if (linsys_is_regtype_ridge (l->lsys)) {	// Ridge
			/*** c -= scale2 * lambda2 * E' * E * beta ***/
			daxpy_ (LINSYS_CINTP (p), &alpha, l->beta, &ione, l->c, &ione);
		} else {
			/*** c -= scale2 * lambda2 * J' * J * beta ***/
			size_t			pj = linsys_get_pj (l->lsys);
			const double	*jr = linsys_get_penalty (l->lsys);
			double			*jb = (double *) malloc (pj * sizeof (double));
			double			*jtjb = (double *) malloc (p * sizeof (double));

			// J * beta
			dgemv_ ("N", LINSYS_CINTP (pj), LINSYS_CINTP (p), &done, jr, LINSYS_CINTP (pj), l->beta, &ione, &dzero, jb, &ione);
			// J' * (J * beta)
			dgemv_ ("T", LINSYS_CINTP (pj), LINSYS_CINTP (p), &done, jr, LINSYS_CINTP (pj), jb, &ione, &dzero, jtjb, &ione);
			free (jb);
			// c -= scale2 * lambda2 * J' * J * beta
			daxpy_ (LINSYS_CINTP (p), &alpha, jtjb, &ione, l->c, &ione);
			free (jtjb);
		}
	}

	{
		int		maxidx = idamax_ (LINSYS_CINTP (p), l->c, &ione) - 1;	// differ of fortran and C
		if (l->sizeA == 0) {
			l->oper.index_of_A = 0;
			l->oper.column_of_X = maxidx;
		}
		l->sup_c = fabs (l->c[maxidx]);
	}

	return;
}

/* Update beta and mu. beta += stepsize * w, mu += stepsze * u */
static void
update_solutions (larsen *l)
{
	size_t		n = linsys_get_n (l->lsys);
	size_t		p = linsys_get_p (l->lsys);
	double		stepsize = (!l->is_interped) ? l->stepsize : l->stepsize_intr;
	double		*beta = (double *) malloc (p * sizeof (double));
	double		*mu = (double *) malloc (n * sizeof (double));

	/*
	 *  in the case of l->interp == true, i.e.,
	 *  to calculate interpolated solution,
	 *  beta_intr = beta_prev + stepsize_intr * w,
	 *  mu_intr = mu_prev + stepsize_intr * u
	 */
	if (!l->is_interped) {
		dcopy_ (LINSYS_CINTP (p), l->beta, &ione, beta, &ione);
		dcopy_ (LINSYS_CINTP (n), l->mu, &ione, mu, &ione);
	} else {
		dcopy_ (LINSYS_CINTP (p), l->beta_prev, &ione, beta, &ione);
		dcopy_ (LINSYS_CINTP (n), l->mu_prev, &ione, mu, &ione);
	}
	larsen_axapy (l, stepsize, l->w, beta);			// beta(A) += stepsize * w(A)
	daxpy_ (LINSYS_CINTP (n), &stepsize, l->u, &ione, mu, &ione);	// mu += stepsize * u

	if (!l->is_interped) {
		dcopy_ (LINSYS_CINTP (p), l->beta, &ione, l->beta_prev, &ione);
		dcopy_ (LINSYS_CINTP (n), l->mu, &ione, l->mu_prev, &ione);
		dcopy_ (LINSYS_CINTP (p), beta, &ione, l->beta, &ione);
		dcopy_ (LINSYS_CINTP (n), mu, &ione, l->mu, &ione);
	} else {
		dcopy_ (LINSYS_CINTP (p), beta, &ione, l->beta_intr, &ione);
		dcopy_ (LINSYS_CINTP (n), mu, &ione, l->mu_intr, &ione);
	}
	free (beta);
	free (mu);

	return;
}

static void
update_stop_loop_flag (larsen *l)
{
	size_t	n = linsys_get_n (l->lsys);
	size_t	p = linsys_get_p (l->lsys);
	int		size = l->sizeA;
	int		m = (linsys_is_regtype_lasso (l->lsys)) ? (int) LINSYS_MIN (n - 1, p) : (int) p;
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
	l->is_interped = false;
	l->stop_loop = true;

	update_correlations (l);

	if (!update_activeset (l)) return false;

	if (!update_equiangular (l)) return false;

	if (!update_stepsize (l)) return false;

	update_solutions (l);

	update_stop_loop_flag (l);

	return true;
}

/* Interpolation
 * In the case of l->lambda1 < | beta | after larsen_regression_step (),
 * the solution corresponding to a designed lambda1 is obtained by the
 * following interpolation:
 * beta_intr = beta_prev + l->stepsize_intr * w
 * mu_intr = mu_prev + l->stepsize_intr * u
 * where l->stepsize_intr = l->absA * (l->lambda1 - | beta_prev |)
 */
bool
larsen_interpolate (larsen *l)
{
	size_t	p = linsys_get_p (l->lsys);
	double	nrm1_prev = dasum_ (LINSYS_CINTP (p), l->beta_prev, &ione);
	double	nrm1 = dasum_ (LINSYS_CINTP (p), l->beta, &ione);
	double	lambda1 = larsen_get_lambda1 (l, true);

	l->is_interped = false;
	if (nrm1_prev <= lambda1 && lambda1 < nrm1) {
		l->is_interped = true;
		l->stepsize_intr = l->absA * (lambda1 - nrm1_prev);
		update_solutions (l);
	}
	return l->is_interped;
}
