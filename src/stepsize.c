/*
 * stepsize.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdlib.h>
#include <math.h>
#include <larsen.h>

#include "linsys_private.h"

/* activeset.c */
extern int		*complementA (larsen *l);

extern void	larsen_awpy (larsen *l, double alpha, double *w, double *y);
extern double	*larsen_xa_dot_ya (larsen *l, const size_t n, double alpha, const double *x, const double *ya);
extern double	*larsen_xa_transpose_dot_y (larsen *l, const size_t n, const double alpha, const double *x, const double *y);

/* \hat{gamma} */
static void
calc_gamma_hat (larsen *l, int *index, int *column, double *val)
{
	int			minplus_idx = -1;
	double		minplus = LINSYS_POSINF;
	int			*Ac = complementA (l);

	if (l->sizeA == l->sys->p) {
		minplus = l->sup_c / l->absA;
	} else if (l->sys->p > l->sizeA) {
		int			i;		double		*a = (double *) malloc (l->sys->p * sizeof (double));
		/* a = scale * X' * u */
		dgemv_ ("T", LINSYS_CINTP (l->sys->n), LINSYS_CINTP (l->sys->p), &l->scale, l->sys->x, LINSYS_CINTP (l->sys->n), l->u, &ione, &dzero, a, &ione);
		if (!l->sys->pen) {	// elastic net
			/* For elastic net, a(A) must be a(A) += l->lambda2 * l->scale^2 * w.
			 * But for the estimation of gamma_hat, a(A) are not referred. */
			/* do nothing */
		} else {
			/* a(A) += l->lambda2 * l->scale^2 * JA' * JA * w */
			size_t	p1 = l->sys->pen->p1;
			double	alpha = l->lambda2 * l->scale2;
			/* jw = JA * w */
			double	*jw = larsen_xa_dot_ya (l, p1, 1., l->sys->pen->r, l->w);
			/* jtjw = JA' * JA * w */
			double	*jtjw = larsen_xa_transpose_dot_y (l, p1, 1., l->sys->pen->r, jw);
			free (jw);
			/* a(A) += alpha * JA' * JA * w */
			larsen_awpy (l, alpha, jtjw, a);
			free (jtjw);
		}
		for (i = 0; i < l->sys->p - l->sizeA; i++) {
			int		j = Ac[i];
			double	cj = l->c[j];
			double	aj = a[j];
			double	e0, e1, min;
			e0 = (l->sup_c - cj) / (l->absA - aj);
			e1 = (l->sup_c + cj) / (l->absA + aj);

			if (e0 <= 0.) e0 = LINSYS_POSINF;
			if (e1 <= 0.) e1 = LINSYS_POSINF;
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
	double	minplus = LINSYS_POSINF;

	if (l->sizeA > 0) {
		int		i;
		for (i = 0; i < l->sizeA; i++) {
			int		j = l->A[i];
			double	e = - l->beta[j] / l->w[i];
			if (e <= 0.) e = LINSYS_POSINF;
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
check_stepsize (const double stepsize)
{
	return (0 < stepsize && stepsize != LINSYS_POSINF);
}

/* Update stepsize and activeset operation l->oper.
 * If gamma_hat > gamma_tilde, the activeset operation on the next
 * step is add a new variable to the active set, else remove a
 * variable from the active set
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
