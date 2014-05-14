/*
 * larsen.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdlib.h>
#include <math.h>
#include <larsen.h>

#define LARSEN_MIN(a, b)		(((a) > (b)) ? (b) : (a))

/* active_set.c */
extern bool	update_activeset (larsen *l);

/* stepsize.c */
extern bool	update_stepsize (larsen *l);

/* equiangular.c */
extern bool	update_equiangular (larsen *l);

/* y(A) = alpha * w(A) + y(A) */
static void
larsen_awpy (larsen *l, double alpha, double *w, double *y)
{
	int		i;
	for (i = 0; i < l->sizeA; i++) {
		int		j = l->A[i];
		y[j] += alpha * w[i];
	}
	return;
}

/* correlation between residual (y - X * beta) and predictors X' */
static void
update_correlations (larsen *l)
{
	double	*r = (double *) malloc (l->n * sizeof (double));
	cblas_dcopy (l->n, l->y, 1, r, 1);
	cblas_daxpy (l->n, -1., l->mu, 1, r, 1);	// r = - mu + y

	/*
	 *  c = Z' * (b - Z * beta),  b = [y; 0], Z = scale * [X; sqrt(lambda2) * E]
	 *  if lambda2 > 0 (scale != 1),
	 *  c = scale * (X' * (y - mu) - scale * lambda2 * beta)
	 */
	if (l->c) free (l->c);
	l->c = (double *) malloc (l->p * sizeof (double));
	cblas_dgemv (CblasColMajor, CblasTrans, l->n, l->p, l->scale, l->x, l->n, r, 1, 0., l->c, 1);
	if (l->is_elnet) cblas_daxpy (l->p, - l->lambda2 * l->scale2, l->beta, 1, l->c, 1);
	free (r);

	{
		int		maxidx = cblas_idamax (l->p, l->c, 1);
		if (l->sizeA == 0) {
			l->oper.index_of_A = 0;
			l->oper.column_of_X = maxidx;
		}
		l->sup_c = fabs (l->c[maxidx]);
	}

	return;
}

/* update beta and mu. beta += stepsize * w, mu += stepsze * u */
static void
update_solutions (larsen *l)
{
	double		stepsize = (!l->interp) ? l->stepsize : l->stepsize_intr;
	double		*beta = (double *) malloc (l->p * sizeof (double));
	double		*mu = (double *) malloc (l->n * sizeof (double));

	/*
	 *  in the case of l->interp == true, i.e.,
	 *  to calculate interpolated solution,
	 *  beta_intr = beta_prev + stepsize_intr * w,
	 *  mu_intr = mu_prev + stepsize_intr * u
	 */
	if (!l->interp) {
		cblas_dcopy (l->p, l->beta, 1, beta, 1);
		cblas_dcopy (l->n, l->mu, 1, mu, 1);
	} else {
		cblas_dcopy (l->p, l->beta_prev, 1, beta, 1);
		cblas_dcopy (l->n, l->mu_prev, 1, mu, 1);
	}
	larsen_awpy (l, stepsize, l->w, beta);			// beta(A) += stepsize * w(A)
	cblas_daxpy (l->n, stepsize, l->u, 1, mu, 1);	// mu += stepsize * u

	if (!l->interp) {
		cblas_dcopy (l->p, l->beta, 1, l->beta_prev, 1);
		cblas_dcopy (l->n, l->mu, 1, l->mu_prev, 1);
		cblas_dcopy (l->p, beta, 1, l->beta, 1);
		cblas_dcopy (l->n, mu, 1, l->mu, 1);
	} else {
		cblas_dcopy (l->p, beta, 1, l->beta_intr, 1);
		cblas_dcopy (l->n, mu, 1, l->mu_intr, 1);
	}
	free (beta);
	free (mu);

	return;
}

static void
update_stop_loop_flag (larsen *l)
{
	int		size = l->sizeA;
	int		n = (l->is_elnet) ? l->p : LARSEN_MIN (l->n - 1, l->p);
	if (l->oper.action == ACTIVESET_ACTION_DROP) size--;
	l->stop_loop = (size >= n) ? true : false;
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
	l->interp = false;
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
	double	lambda1 = (l->is_elnet) ? l->scale * l->lambda1 : l->lambda1;	// scale * lambda1
	double	nrm1_prev = cblas_dasum (l->p, l->beta_prev, 1);
	double	nrm1 = cblas_dasum (l->p, l->beta, 1);
	l->interp = false;
	if (nrm1_prev <= lambda1 && lambda1 < nrm1) {
		l->interp = true;
		l->stepsize_intr = l->absA * (lambda1 - nrm1_prev);
		update_solutions (l);
	}
	return l->interp;
}
