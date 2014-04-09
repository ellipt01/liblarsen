/*
 * larsen.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <larsen.h>

#define MIN(a, b)		((a > b) ? b : a)

/* active_set.c */
extern bool	activeset_add (larsen *l, int j);
extern bool	activeset_remove (larsen *l, int j);

/* stepsize.c */
extern bool	update_stepsize (larsen *l);

/* equiangular.c */
extern bool	update_equiangular (larsen *l);

/* y(A) = alpha * w(A) + y(A) */
void
larsen_awpy (larsen *l, double alpha, c_vector *w, c_vector *y)
{
	int				i;
	for (i = 0; i < l->A->size; i++) {
		int		j = c_vector_int_get (l->A, i);
		double	wi = c_vector_get (w, i);
		double	yj = c_vector_get (y, j);
		c_vector_set (y, j, yj + alpha * wi);
	}
	return;
}

/* correlation between residual (y - X * beta) and predictors X' */
static void
update_correlations (larsen *l)
{
	c_vector	*r = c_vector_alloc (l->y->size);
	c_vector_memcpy (r, l->y);
	c_vector_sub (r, l->mu);

	/*
	 *  c = Z' * (b - Z * beta),  b = [y; 0], Z = scale * [X; sqrt(lambda2) * E]
	 *  if lambda2 > 0 (scale != 1),
	 *  c = scale * (X' * (y - mu) - scale * lambda2 * beta)
	 */
	if (l->c) c_vector_free (l->c);
	l->c = c_matrix_transpose_dot_vector (1., l->x, r, 0.);
	if (l->is_elnet) {
		c_vector_axpy (- l->lambda2 * l->scale, l->beta, l->c);
		c_vector_scale (l->c, l->scale);
	}
	c_vector_free (r);

	{
		int		maxidx = c_vector_amax (l->c);
		if (l->A->size == 0) l->oper.column = maxidx;
		l->sup_c = fabs (c_vector_get (l->c, maxidx));
	}

	return;
}

/* update beta and mu. beta += stepsize * w, mu += stepsze * u */
static void
update_solutions (larsen *l)
{
	double		stepsize = (!l->interp) ? l->stepsize : l->stepsize_intr;
	c_vector	*beta = c_vector_alloc (l->beta->size);
	c_vector	*mu = c_vector_alloc (l->mu->size);

	/*
	 *  in the case of l->interp == true, i.e.,
	 *  to calcurate interpolated solution,
	 *  beta_intr = beta_prev + stepsize_intr * w,
	 *  mu_intr = mu_prev + stepsize_intr * u
	 */
	if (!l->interp) {
		c_vector_memcpy (beta, l->beta);
		c_vector_memcpy (mu, l->mu);
	} else {
		c_vector_memcpy (beta, l->beta_prev);
		c_vector_memcpy (mu, l->mu_prev);
	}
	larsen_awpy (l, stepsize, l->w, beta);
	c_vector_axpy (stepsize, l->u, mu);

	if (!l->interp) {
		c_vector_memcpy (l->beta_prev, l->beta);
		c_vector_memcpy (l->mu_prev, l->mu);
		c_vector_memcpy (l->beta, beta);
		c_vector_memcpy (l->mu, mu);
	} else {
		c_vector_memcpy (l->beta_intr, beta);
		c_vector_memcpy (l->mu_intr, mu);
	}
	c_vector_free (beta);
	c_vector_free (mu);

	return;
}

static void
update_stop_loop_flag (larsen *l)
{
	int		size = l->A->size;
	int		n = (l->is_elnet) ? l->x->size2 : MIN (l->x->size1 - 1, l->x->size2);
	if (l->oper.action == ACTIVESET_ACTION_DROP) size--;
	l->stop_loop = (size >= n) ? true : false;
	if (!l->stop_loop) l->stop_loop = (l->oper.column == -1);
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

	if (l->oper.action == ACTIVESET_ACTION_ADD) activeset_add (l, l->oper.column);
	else if (l->oper.action == ACTIVESET_ACTION_DROP) activeset_remove (l, l->oper.column);
	else return false;

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
	double	nrm1_prev = c_vector_asum (l->beta_prev);
	double	nrm1 = c_vector_asum (l->beta);
	l->interp = false;
	if (nrm1_prev <= lambda1 && lambda1 < nrm1) {
		l->interp = true;
		l->stepsize_intr = l->absA * (lambda1 - nrm1_prev);
		update_solutions (l);
	}
	return l->interp;
}

/* return copy of elastic net solution: beta_elnet = scale * l->beta */
c_vector *
larsen_get_beta (larsen *l)
{
	c_vector	*beta = c_vector_alloc (l->beta->size);
	if (!l->interp) c_vector_memcpy (beta, l->beta);
	else c_vector_memcpy (beta, l->beta_intr);
	if (l->is_elnet) c_vector_scale (beta, 1. / l->scale);
	return beta;
}

/* return copy of elastic net solution: mu_elnet = scale^2 * mu_navie */
c_vector *
larsen_get_mu (larsen *l)
{
	c_vector	*mu = c_vector_alloc (l->mu->size);
	if (!l->interp) c_vector_memcpy (mu, l->mu);
	else c_vector_memcpy (mu, l->mu_intr);
	if (l->is_elnet) c_vector_scale (mu, pow (1. / l->scale, 2.));
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
