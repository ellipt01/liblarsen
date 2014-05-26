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

/* evaluate min = MIN_+ (a, b)
 * if min <= 0, return false */
static bool
eval_minplus (const double a, const double b, double *min)
{
	if (a <= 0.) *min = b;
	else if (b <= 0.) *min = a;
	else *min = (a <= b) ? a : b;

	return (*min > 0.);
}

/* gamma_hat */
static void
calc_gamma_hat (larsen *l, int *index, int *column, double *val)
{
	size_t		n = l->lreg->n;
	size_t		p = l->lreg->p;
	int			minplus_idx = -1;
	double		minplus = -1.;
	int			*Ac = complementA (l);

	if (l->sizeA == p) {
		minplus = l->sup_c / l->absA;
	} else if (p > l->sizeA) {
		int				i;
		double			scale = l->lreg->scale;
		const double	*x = l->lreg->x;
		double			*a = (double *) malloc (p * sizeof (double));

		/*** a = scale * X' * u ***/
		dgemv_ ("T", LINREG_CINTP (n), LINREG_CINTP (p), &scale, x, LINREG_CINTP (n), l->u, &ione, &dzero, a, &ione);

		/*** a += l->lambda2 * l->scale^2 * J' * J(:,A) * w ***/
		if (!larsen_is_regtype_lasso (l)) {

			if (larsen_is_regtype_ridge (l)) {	// Ridge
				/* For elastic net,
				 *   a += z
				 *   z = l->lambda2 * l->scale^2 * E' * E(:,A) * w,
				 * and z(A) != 0, z(Ac) = 0.
				 * But for the estimation of gamma_hat, a(A), and z(A) are not referred. */

				/* do nothing */
			} else {
				/*** a += l->lambda2 * l->scale^2 * J' * J(:,A) * w ***/
				size_t			pj = l->lreg->pen->pj;
				const double	*jr = l->lreg->pen->r;
				double			alpha = l->lreg->lambda2 * l->lreg->scale2;
				// J(:,A) * w
				double			*jw = larsen_xa_dot_ya (l, pj, 1., jr, l->w);
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
			double	eminus, eplus, min;
			eminus = (l->sup_c - cj) / (l->absA - aj);
			eplus = (l->sup_c + cj) / (l->absA + aj);

			if (!eval_minplus (eminus, eplus, &min)) continue;
			if (minplus < 0. || min < minplus) {
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

/* gamma_tilde */
static void
calc_gamma_tilde (larsen *l, int *index, int *column, double *val)
{
	int		minplus_idx = -1;
	double	minplus = -1.;

	if (l->sizeA > 0) {
		int		i;
		for (i = 0; i < l->sizeA; i++) {
			int		j = l->A[i];
			double	e = - l->beta[j] / l->w[i];
			if (e <= 0.) continue;
			if (minplus < 0. || e < minplus) {
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
	return (stepsize > 0.);
}

enum {
	GAMMA_HAT_SELECTED	= 1,
	GAMMA_TILDE_SELECTED	= 2
};

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

	int		gamma_status = 0;	// GAMMA_HAT_SELECTED or GAMMA_TILDE_SELECTED

	l->stepsize = -1.;
	l->oper.action = ACTIVESET_ACTION_NONE;
	l->oper.index_of_A = -1;
	l->oper.column_of_X = -1;

	calc_gamma_hat (l, &gamma_hat_idx, &gamma_hat_col, &gamma_hat);
	calc_gamma_tilde (l, &gamma_tilde_idx, &gamma_tilde_col, &gamma_tilde);

	if (gamma_hat <= 0. && gamma_tilde <= 0.) return false;
	if (gamma_hat <= 0.) gamma_status = GAMMA_TILDE_SELECTED;
	else if (gamma_tilde <= 0.) gamma_status = GAMMA_HAT_SELECTED;
	else {
		/* if gamma_tilde == gamma_hat, gamma = gamma_tilde (remove a predictor from active set) */
		gamma_status = (gamma_tilde <= gamma_hat) ? GAMMA_TILDE_SELECTED : GAMMA_HAT_SELECTED;
	}

	if (gamma_status == GAMMA_HAT_SELECTED) {
		l->oper.action = ACTIVESET_ACTION_ADD;
		l->stepsize = gamma_hat;
		if (gamma_hat_idx >= 0) l->oper.index_of_A = gamma_hat_idx;
		if (gamma_hat_col >= 0) l->oper.column_of_X = gamma_hat_col;
	} else if (gamma_status == GAMMA_TILDE_SELECTED) {	// gamma_tilde <= gamma_hat
		l->oper.action = ACTIVESET_ACTION_DROP;
		l->stepsize = gamma_tilde;
		if (gamma_tilde_idx >= 0) l->oper.index_of_A = gamma_tilde_idx;
		if (gamma_tilde_col >= 0) l->oper.column_of_X = gamma_tilde_col;
	}

	return check_stepsize (l->stepsize);
}
