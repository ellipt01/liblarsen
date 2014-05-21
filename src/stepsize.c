/*
 * stepsize.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdlib.h>
#include <math.h>
#include <larsen.h>

#include "larsen_private.h"

/* \hat{gamma} */
static void
calc_gamma_hat (larsen *l, int *index, int *column, double *val)
{
	size_t		n = linreg_get_n (l->lreg);
	size_t		p = linreg_get_p (l->lreg);
	int			minplus_idx = -1;
	double		minplus = LINREG_POSINF;	// (+1. / 0.)
	int			*Ac = complementA (l);

	if (l->sizeA == p) {
		minplus = l->sup_c / l->absA;
	} else if (p > l->sizeA) {
		int				i;
		double			scale = linreg_get_scale (l->lreg);
		const double	*x = linreg_get_x (l->lreg);
		double			*a = (double *) malloc (p * sizeof (double));

		/*** a = scale * X' * u ***/
		dgemv_ ("T", LINREG_CINTP (n), LINREG_CINTP (p), &scale, x, LINREG_CINTP (n), l->u, &ione, &dzero, a, &ione);

		/*** a += l->lambda2 * l->scale^2 * J' * J(:,A) * w ***/
		if (!linreg_is_regtype_lasso (l->lreg)) {

			if (linreg_is_regtype_ridge (l->lreg)) {	// Ridge
				/* For elastic net,
				 *   a += z
				 *   z = l->lambda2 * l->scale^2 * E' * E(:,A) * w,
				 * and z(A) != 0, z(Ac) = 0.
				 * But for the estimation of gamma_hat, a(A), and z(A) are not referred. */

				/* do nothing */
			} else {
				/*** a += l->lambda2 * l->scale^2 * J' * J(:,A) * w ***/
				size_t			pj = linreg_get_pj (l->lreg);
				const double	*jr = linreg_get_penalty (l->lreg);
				double			alpha = linreg_get_lambda2 (l->lreg) * linreg_get_scale2 (l->lreg);
				double			*jw = larsen_xa_dot_ya (l, pj, 1., jr, l->w);	// J(:,A) * w

				// J' * (J(:,A) * w)
				double			*jtjw = (double *) malloc (p * sizeof (double));
				dgemv_ ("T", LINREG_CINTP (pj), LINREG_CINTP (p), &done, jr, LINREG_CINTP (pj), jw, &ione, &dzero, jtjw, &ione);
				free (jw);
				/* a += alpha * J' * J(:,A) * w */
				daxpy_ (LINREG_CINTP (p), &alpha, jtjw, &ione, a, &ione);
				free (jtjw);

			}
		}

		for (i = 0; i < p - l->sizeA; i++) {
			int		j = Ac[i];
			double	cj = l->c[j];
			double	aj = a[j];
			double	e0, e1, min;
			e0 = (l->sup_c - cj) / (l->absA - aj);
			e1 = (l->sup_c + cj) / (l->absA + aj);

			if (e0 <= 0.) e0 = LINREG_POSINF;
			if (e1 <= 0.) e1 = LINREG_POSINF;
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
	double	minplus = LINREG_POSINF;

	if (l->sizeA > 0) {
		int		i;
		for (i = 0; i < l->sizeA; i++) {
			int		j = l->A[i];
			double	e = - l->beta[j] / l->w[i];
			if (e <= 0.) e = LINREG_POSINF;
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
	return (0 < stepsize && stepsize != LINREG_POSINF);
}
int count = 0;
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
