/*
 * stepsize.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdlib.h>
#include <larsen.h>

#define posinf()	(1. / 0.)	// + infinity

/* activeset.c */
extern int		*complementA (larsen *l);

/* \hat{gamma} */
static void
calc_gamma_hat (larsen *l, int *index, int *column, double *val)
{
	int			minplus_idx = -1;
	double		minplus = posinf ();
	int			*Ac = complementA (l);

	if (l->sizeA == l->p) {
		minplus = l->sup_c / l->absA;
	} else if (l->p > l->sizeA) {
		int			i;
		double		*a = (double *) malloc (l->p * sizeof (double));
		cblas_dgemv (CblasColMajor, CblasTrans, l->n, l->p, l->scale, l->x, l->n, l->u, 1, 0., a, 1);
		for (i = 0; i < l->p - l->sizeA; i++) {
			int		j = Ac[i];
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
	*index = (minplus_idx >= 0) ? l->sizeA : -1;
	*column = (minplus_idx >= 0) ? Ac[minplus_idx] : -1;
	*val = minplus;
	free (Ac);
	return;
}

/* \tilde{gamma} */
static void
calc_gamma_tilde (larsen *l, int *index, int *column, double *val)
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
	*index = (minplus_idx >= 0) ? minplus_idx : -1;
	*column = (minplus_idx >= 0) ? l->A[minplus_idx] : -1;
	*val = minplus;
	return;
}

static bool
check_stepsize (double stepsize)
{
	return (stepsize != posinf ());
}

/*
 * if gamma_hat > gamma_tilde, add a new variable to the active set
 * else remove a variable from the active set
 */
bool
update_stepsize (larsen *l)
{
	int		gamma_hat_idx;
	int		gamma_hat_col;
	double	gamma_hat;
	int		gamma_tilde_idx;
	int		gamma_tilde_col;
	double	gamma_tilde;

	l->oper.action = ACTIVESET_ACTION_NONE;

	calc_gamma_hat (l, &gamma_hat_idx, &gamma_hat_col, &gamma_hat);
	calc_gamma_tilde (l, &gamma_tilde_idx, &gamma_tilde_col, &gamma_tilde);

	l->oper.index_of_A = -1;
	l->oper.column_of_X = -1;
	if (gamma_hat < gamma_tilde) {
		l->oper.action = ACTIVESET_ACTION_ADD;
		l->stepsize = gamma_hat;
		if (gamma_hat_idx >= 0) l->oper.index_of_A = gamma_hat_idx;
		if (gamma_hat_col >= 0) l->oper.column_of_X = gamma_hat_col;
	} else {	// gamma_tilde <= gamma_hat
		l->oper.action = ACTIVESET_ACTION_DROP;
		l->stepsize = gamma_tilde;
		if (gamma_tilde_idx >= 0) l->oper.index_of_A = gamma_tilde_idx;
		if (gamma_tilde_col >= 0) l->oper.column_of_X = gamma_tilde_col;
	}

	return check_stepsize (l->stepsize);
}
