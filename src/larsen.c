/*
 * larsen.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdlib.h>
#include <math.h>
#include <larsen.h>

#define MIN(a, b)		(((a) > (b)) ? (b) : (a))

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
	int		n = (l->is_elnet) ? l->p : MIN (l->n - 1, l->p);
	if (l->oper.action == ACTIVESET_ACTION_DROP) size--;
	l->stop_loop = (size >= n) ? true : false;
	if (!l->stop_loop) l->stop_loop = (l->oper.column_of_X == -1);
	return;
}

/* Progress one step of the LARS-EN algorithm
 * add / remove one variable (assigned by l->oper.column)
 * to / from the active set.
 *
 * 1. Update correlation between residualrelevants (y - mu) and variables (X').
 *    In the case of A = {} (first iteration or restart), a variable
 *    which has largest correlation is selected (l->column is set to its index in X).
 *
 * 2. Add / remove one variable assigned by l->oper.column to / from the active set.
 *
 * 3. Update equi-angular vector and its relevant.
 *
 * 4. Update step size and l->oper (which specify the next operation to the active set).
 *    l->oper.action is updated according to whether gamma_hat < gamma_tilde or not.
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
	double	lambda1 = larsen_get_lambda1 (l, true);
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

/* return copy of elastic net solution: beta_elnet = scale * l->beta */
double *
larsen_get_beta (larsen *l)
{
	double	*beta = (double *) malloc (l->p * sizeof (double));
	if (!l->interp) cblas_dcopy (l->p, l->beta, 1, beta, 1);
	else cblas_dcopy (l->p, l->beta_intr, 1, beta, 1);
	if (l->is_elnet) cblas_dscal (l->p, 1. / l->scale, beta, 1);
	return beta;
}

/* return copy of elastic net solution: mu_elnet = scale^2 * mu_navie */
double *
larsen_get_mu (larsen *l)
{
	double	*mu = (double *) malloc (l->n * sizeof (double));
	if (!l->interp) cblas_dcopy (l->n, l->mu, 1, mu, 1);
	else cblas_dcopy (l->n, l->mu_intr, 1, mu, 1);
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

/* return l->lambda1, if (scaled && l->is_elnet) return l->lambda1 * l->scale */
double
larsen_get_lambda1 (larsen *l, const bool scaled)
{
	if (!scaled) return l->lambda1;
	return (l->is_elnet) ? l->lambda1 * l->scale : l->lambda1;
}
