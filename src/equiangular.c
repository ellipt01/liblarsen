/*
 * equiangular.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <larsen.h>

#include "larsen_private.h"

/* s = sign(c) */
static double *
update_sign (larsen *l)
{
	int			i;
	double		*s = (double *) malloc (l->sizeA * sizeof (double));
	for (i = 0; i < l->sizeA; i++) {
		int		j = l->A[i];
		s[i] = (l->c[j] >= 0.) ? 1. : -1.;
	}
	return s;
}

/* check result of Cholesky insert or delete */
static bool
check_info (const char *func_name, const int info)
{
	bool	is_valid_info = false;

	if (info == 0) is_valid_info = true;
	else if (info < 0) fprintf (stderr, "WARNING: %s : the %d-th argument had an illegal value.\n", func_name, - info - 1);
	else if (info > 0) {
		fprintf (stderr, "WARNING: %s : the leading minor of order %d is not positive definite.\n", func_name, info);
		fprintf (stderr, "         factorization could not be completed.\n");
	}
	return is_valid_info;
}

/* do cholinsert / choldelete of specified index */
static int
update_chol (larsen *l)
{
	int		info = 0;
	int		index = l->oper.index_of_A;

	if (l->oper.action == ACTIVESET_ACTION_ADD) {
		/*** insert a predictor ***/
		int				j = l->oper.column_of_X;
		double			*t = (double *) malloc (l->sizeA * sizeof (double));
		size_t			n = l->lreg->n;
		const double	*x = l->lreg->x;
		const double	*xj = x + LINREG_INDEX_OF_MATRIX (0, j, n);
		double			scale2 = l->lreg->scale2;

		/*** t = scale^2 * X(:,A)' * X(:,j) ***/
		t = larsen_xa_transpose_dot_y (l, n, scale2, x, xj);

		if (!larsen_is_regtype_lasso (l)) {
			double		lambda2 = l->lreg->lambda2;
			double		alpha = lambda2 * scale2;
			/* In the case of lambda2 > 0, (now, A is already updated and j \in A)
			 *
			 *   t = Z(:,A)' * Z(:,j)
			 *     = scale^2 * [X(:,A)', sqrt(lambda2) * J(:,A)'] * [X(:,j); sqrt(lambda2) * J(:,j)]
			 *     = scale^2 * [X(:,A)' * X(:,j) + lambda2 * J(:,A)' * J(:,j)]
			 * So,
			 *   t += scale^2 * lambda2 * J(:,A)' * J(:,j)
			 */
			if (larsen_is_regtype_ridge (l)) {	// Ridge
				/* In this case, J = E, so
				 *   t += scale^2 * lambda2 * E(:,A)' * E(:,j)
				 *      = scale^2 * [E(:,A-1)' * E(:,j) (= 0); E(:,j)' * E(:,j) (= 1)]
				 * so,
				 *   t[l->oper.index_of_A] += scale^2 * lambda2
				 * where, because X(:,j)' * X(:,j) = 1,
				 * t[l->oper.index_of_A] = scale^2 + lambda2 * scale^2 = 1,
				 * but in the case of adaptive elastic net, the above != 1 */
				t[index] += alpha;
			} else {
				/*** t += scale^2 * lambda2 * J(:,A)' * J(:,j) ***/
				size_t			pj = l->lreg->pen->pj;
				const double	*jr = l->lreg->pen->r;	// J
				const double	*jrj = jr + LINREG_INDEX_OF_MATRIX (0, j, pj);	// J(:,j)
				// J(:,A)' * J(:,j)
				double			*jtj = larsen_xa_transpose_dot_y (l, pj, 1., jr, jrj);
				/* t += scale^2 * lambda2 * J(:,A)' * J(:,j) */
				daxpy_ (LINREG_CINTP (l->sizeA), &alpha, jtj, &ione, t, &ione);
				free (jtj);
			}
		}
		/*** insert t ***/
		info = larsen_linalg_cholesky_insert (l->sizeA - 1, &l->chol, index, t);
		free (t);

	} else if (l->oper.action == ACTIVESET_ACTION_DROP) {
		/*** delete a predictor ***/
		larsen_linalg_cholesky_delete (l->sizeA + 1, &l->chol, index);
	}

	return info;
}

/* update equiangular vector (l->u) and its relevant (l->w and l->absA)
 * for current active set using cholinsert / delete routine */
static bool
update_equiangular_larsen_cholesky (larsen *l)
{
	int			info;
	double		*s;

	if (l->sizeA <= 0) return false;

	/* s = sign (c) */
	s = update_sign (l);

	if (l->w) free (l->w);
	l->w = (double *) malloc (l->sizeA * sizeof (double));
	dcopy_ (LINREG_CINTP (l->sizeA), s, &ione, l->w, &ione);

	/* Cholesky update and solve equiangular equation
	 * TODO: create new methods to solve equiangular equation (using QR, SVD etc.)
	 * and switch them in the following line
	 */
	info = update_chol (l);
	if (!check_info ("update_chol", info)) return false;

	info = larsen_linalg_cholesky_svx (l->sizeA, l->chol, l->sizeA, l->w);
	if (!check_info ("cholesky_svx", info)) return false;
	/* end */

	l->absA = 1. / sqrt (ddot_ (LINREG_CINTP (l->sizeA), s, &ione, l->w, &ione));
	free (s);

	/* w = (ZA' * ZA)^-1 * s(A) * absA */
	dscal_ (LINREG_CINTP (l->sizeA), &l->absA, l->w, &ione);

	/* In exactly
	 *
	 *     u = Z(:,A) * w = scale * [X(:,A) * w; sqrt(lambda2) * J(:,A) * w],
	 *     dim(u) = n + pJ (size of nonzero J(:,A) * w is |A|).
	 *
	 * But in this program, to save memory, only first n part of u is stored
	 *
	 *     u = scale * XA * w,
	 *     dim(u) = n,
	 *
	 * because u(n+1:n+pJ) (= scale * sqrt(lambda2) * J * w ) is not appear
	 * in the positive on the LARS-EN algorithm.
	 */
	if (l->u) free (l->u);
	{
		size_t			n = l->lreg->n;
		double			scale = l->lreg->scale;
		const double	*x = l->lreg->x;
		l->u = larsen_xa_dot_ya (l, n, scale, x, l->w);
	}

	return true;
}

/* call equiangular vector updater
 * if another routine is implemented, switch them in here */
bool
update_equiangular (larsen *l)
{
	return update_equiangular_larsen_cholesky (l);
}
