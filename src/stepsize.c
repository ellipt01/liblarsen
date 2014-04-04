/*
 * stepsize.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <larsen.h>

extern void	larsen_awpy (larsen *l, double alpha, cl_vector *w, cl_vector *y);

static double
posinf (void)
{
	return + 1. / 0.;
}

/* \hat{gamma} */
static void
calc_gamma_hat (larsen *l, int *index, double *val)
{
	int			minplus_idx = -1;
	double		minplus = posinf ();

	if (l->A->size == l->x->size2) {

		minplus = l->sup_c / l->absA;

	} else if (l->Ac->size > 0) {
		int			i;
		cl_vector	*a = cl_matrix_transpose_dot_vector (l->x, l->u);
		if (l->do_scaling) {
			larsen_awpy (l, l->lambda2 * l->scale, l->w, a);
			cl_vector_scale (a, l->scale);
		}
		for (i = 0; i < l->Ac->size; i++) {
			int		j = cl_vector_int_get (l->Ac, i);
			double	cj = cl_vector_get (l->c, j);
			double	aj = cl_vector_get (a, j);
			double	e0, e1, min;
			e0 = (l->sup_c - cj) / (l->absA - aj);
			e1 = (l->sup_c + cj) / (l->absA + aj);

			if (e0 <= 0.) e0 = posinf ();
			if (e1 <= 0.) e1 = posinf ();
			min = (e0 <= e1) ? e0 : e1;

			if (min < minplus) {
				minplus_idx = i;
				minplus = min;
			}
		}
		cl_vector_free (a);
	}
	*index = minplus_idx;
	*val = minplus;
	return;
}

/* \tilde{gamma} */
static void
calc_gamma_tilde (larsen *l, int *index, double *val)
{
	int		minplus_idx = -1;
	double	minplus = posinf ();

	if (l->A->size > 0) {
		int		i;
		for (i = 0; i < l->A->size; i++) {
			int		j = cl_vector_int_get (l->A, i);
			double	betaj = cl_vector_get (l->beta, j);
			double	wi = cl_vector_get (l->w, i);
			double	e = - betaj / wi;
			if (e <= 0.) e = posinf ();
			if (e < minplus) {
				minplus_idx = i;
				minplus = e;
			}
		}
	}
	*index = minplus_idx;
	*val = minplus;
	return;
}

/*
 * if gamma_hat > gamma_tilde, add a new variable to the active set
 * else remove a variable from the active set
 */
bool
update_stepsize (larsen *l)
{
	l->oper.action = ACTIVESET_ACTION_NONE;
	{
		int		gamma_hat_idx;
		double	gamma_hat;
		int		gamma_tilde_idx;
		double	gamma_tilde;

		calc_gamma_hat (l, &gamma_hat_idx, &gamma_hat);
		calc_gamma_tilde (l, &gamma_tilde_idx, &gamma_tilde);

		l->oper.column = -1;
		if (gamma_hat < gamma_tilde) {
			l->oper.action = ACTIVESET_ACTION_ADD;
			l->stepsize = gamma_hat;
			if (gamma_hat_idx >= 0) l->oper.column = cl_vector_int_get (l->Ac, gamma_hat_idx);
		} else {
			l->oper.action = ACTIVESET_ACTION_DROP;
			l->stepsize = gamma_tilde;
			if (gamma_tilde_idx >= 0) l->oper.column = cl_vector_int_get (l->A, gamma_tilde_idx);
		}
	}
	return (l->stepsize != posinf ());
}
