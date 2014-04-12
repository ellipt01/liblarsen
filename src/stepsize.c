/*
 * stepsize.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdlib.h>
#include <larsen.h>

/* larsen.c */
extern void	larsen_awpy (larsen *l, double alpha, double *w, double *y);

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

	if (l->sizeA == l->p) {
		minplus = l->sup_c / l->absA;
	} else if (l->p > l->sizeA) {
		int			i;
		double		*a = (double *) malloc (l->p * sizeof (double));
		cblas_dgemv (CblasColMajor, CblasTrans, l->n, l->p, l->scale, l->x, l->n, l->u, 1, 0., a, 1);
		if (l->is_elnet) larsen_awpy (l, l->lambda2 * l->scale2, l->w, a);
		for (i = 0; i < l->p - l->sizeA; i++) {
			int		j = l->Ac[i];
			double	cj = l->c[j];
			double	aj = a[j];
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
		free (a);
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

	if (l->sizeA > 0) {
		int		i;
		for (i = 0; i < l->sizeA; i++) {
			int		j = l->A[i];
			double	e = - l->beta[j] / l->w[i];
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
			if (gamma_hat_idx >= 0) l->oper.column = l->Ac[gamma_hat_idx];
		} else {
			l->oper.action = ACTIVESET_ACTION_DROP;
			l->stepsize = gamma_tilde;
			if (gamma_tilde_idx >= 0) l->oper.column = l->A[gamma_tilde_idx];
		}
	}
	return (l->stepsize != posinf ());
}
