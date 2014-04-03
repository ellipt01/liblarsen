/*
 * larsen.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include "larsen.h"

/* active_set.c */
extern bool	activeset_add (larsen *l, int j);
extern bool	activeset_remove (larsen *l, int j);

/* step_size.c */
extern bool	update_stepsize (larsen *l);

/* equiangular.c */
extern bool	update_equiangular (larsen *l);

/* y(A) = alpha * w(A) + y(A) */
void
larsen_awpy (larsen *l, double alpha, cl_vector *w, cl_vector *y)
{
	int				i;
	for (i = 0; i < l->A->size; i++) {
		int		j = cl_vector_int_get (l->A, i);
		double	wi = cl_vector_get (w, i);
		double	yj = cl_vector_get (y, j);
		cl_vector_set (y, j, yj + alpha * wi);
	}
	return;
}

/* correlation between residual (y - X * beta) and predictors X' */
static void
update_correlations (larsen *l)
{
	cl_vector	*r = cl_vector_alloc (l->y->size);
	cl_vector_memcpy (r, l->y);
	cl_vector_sub (r, l->mu);

	/*
	 *  c = Z' * (b - Z * beta),  b = [y; 0], Z = scale * [X; sqrt(lambda2) * E]
	 *  if lambda2 > 0 (scale != 1),
	 *  c = scale * (X' * (y - mu) - scale * lambda2 * beta)
	 */
	if (l->c) cl_vector_free (l->c);
	l->c = cl_matrix_transpose_dot_vector (l->x, r);
	if (l->do_scaling) {
		cl_vector_axpy (- l->lambda2 * l->scale, l->beta, l->c);
		cl_vector_scale (l->c, l->scale);
	}
	cl_vector_free (r);

	{
		int		maxidx = cl_vector_amax (l->c);
		if (l->A->size == 0) l->oper.column = maxidx;
		l->sup_c = fabs (cl_vector_get (l->c, maxidx));
	}

	return;
}

/* update beta and mu. beta += stepsize * w, mu += stepsze * u */
static void
update_solutions (larsen *l)
{
	double		stepsize = (!l->interp) ? l->stepsize : l->stepsize_intr;
	cl_vector	*beta = cl_vector_alloc (l->beta->size);
	cl_vector	*mu = cl_vector_alloc (l->mu->size);

	/*
	 *  in the case of l->interp == true, i.e.,
	 *  to calcurate interpolated solution,
	 *  beta_intr = beta_prev + stepsize_intr * w,
	 *  mu_intr = mu_prev + stepsize_intr * u
	 */
	if (!l->interp) {
		cl_vector_memcpy (beta, l->beta);
		cl_vector_memcpy (mu, l->mu);
	} else {
		cl_vector_memcpy (beta, l->beta_prev);
		cl_vector_memcpy (mu, l->mu_prev);
	}
	larsen_awpy (l, stepsize, l->w, beta);
	cl_vector_axpy (stepsize, l->u, mu);

	if (!l->interp) {
		cl_vector_memcpy (l->beta_prev, l->beta);
		cl_vector_memcpy (l->mu_prev, l->mu);
		cl_vector_memcpy (l->beta, beta);
		cl_vector_memcpy (l->mu, mu);
	} else {
		cl_vector_memcpy (l->beta_intr, beta);
		cl_vector_memcpy (l->mu_intr, mu);
	}
	cl_vector_free (beta);
	cl_vector_free (mu);

	return;
}

static void
update_stop_loop_flag (larsen *l)
{
	int		size = l->A->size;
	int		n = (l->do_scaling) ? l->x->size2 : CL_MIN (l->x->size1 - 1, l->x->size2);
	if (l->oper.action == ACTIVESET_ACTION_DROP) size--;
	l->stop_loop = (size >= n) ? true : false;
	if (!l->stop_loop) l->stop_loop = (l->oper.column == -1);
	return;
}

/* progress one step of the regression */
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

double
larsen_get_lambda1 (larsen *l)
{
	return (l->do_scaling) ? l->lambda1 * l->scale : l->lambda1;
}

/* interpolation */
bool
larsen_interpolate (larsen *l)
{
	double	lambda1 = larsen_get_lambda1 (l);
	double	nrm1_prev = cl_vector_asum (l->beta_prev);
	double	nrm1 = cl_vector_asum (l->beta);
	l->interp = false;
	if (nrm1_prev <= lambda1 && lambda1 < nrm1) {
		l->interp = true;
		l->stepsize_intr = l->absA * (lambda1 - nrm1_prev);
		update_solutions (l);
	}
	return l->interp;
}

/* get elastic net solution: beta_navie -> beta_elnet = scale * beta_navie */
cl_vector *
larsen_get_beta (larsen *l)
{
	cl_vector	*beta = cl_vector_alloc (l->beta->size);
	if (!l->interp) cl_vector_memcpy (beta, l->beta);
	else cl_vector_memcpy (beta, l->beta_intr);
	if (l->do_scaling) cl_vector_scale (beta, 1. / l->scale);
	return beta;
}

/* get elastic net solution: mu_navie -> mu_elnet = scale^2 * mu_navie */
cl_vector *
larsen_get_mu (larsen *l)
{
	cl_vector	*mu = cl_vector_alloc (l->mu->size);
	if (!l->interp) cl_vector_memcpy (mu, l->mu);
	else cl_vector_memcpy (mu, l->mu_intr);
	if (l->do_scaling) cl_vector_scale (mu, pow (1. / l->scale, 2.));
	return mu;
}

void
larsen_increment_lambda1 (larsen *l, double dt)
{
	l->lambda1 += dt;
	return;
}

bool
larsen_loop_continue (larsen *l, double stop)
{
	if (stop > 0.) return (l->lambda1 <= stop);
	if (!l->stop_loop) return true;
	return false;
}
