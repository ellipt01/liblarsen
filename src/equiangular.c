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

/* check result of cholesky insert or delete */
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
		size_t			n = linsys_get_n (l->lsys);
		const double	*x = linsys_get_x (l->lsys);
		const double	*xj = x + LINSYS_INDEX_OF_MATRIX (0, j, n);
		double			lambda2 = linsys_get_lambda2 (l->lsys);
		double			scale2 = linsys_get_scale2 (l->lsys);
		double			alpha = lambda2 * scale2;

		/* t = scale^2 * X(:,A)' * X(:,j) */
		t = larsen_xa_transpose_dot_y (l, n, scale2, x, xj);

		if (!linsys_is_regtype_lasso (l->lsys)) {	// lambda2 > 0
			/* Now, A is already updated and j \in A.
			 *
			 * t = Z(:,A)' * Z(:,j)
			 * = scale^2 * [X(:,A)', sqrt(lambda2) * J(:,A)'] * [X(:,j); sqrt(lambda2) * J(:,j)]
			 * = scale^2 * X(:,A)' * X(:,j) + scale^2 * lambda2 * J(:,A) * J(:,j) = 0)
			 *
			 * In the case of elastic net (J = E)
			 * t = scale^2 * [X(:,A)', sqrt(lambda2) * E(:,A)'] * [X(:,j); sqrt(lambda2) * E(:,j)]
			 * = scale^2 * X(:,A-1)' * X(:,j)
			 *   scale^2 * X(:,j)' * X(:,j) + scale^2 * lambda2 * E(:,j) * E(:,j)
			 * so,
			 * t[l->oper.index_of_A] += scale^2 * lambda2
			 * where, because X(:,j)' * X(:,j) = 1,
			 * t[l->oper.index_of_A] = scale^2 + lambda2 * scale^2 = 1,
			 * but in the case of adaptive elastic net, the above != 1 */
			if (linsys_is_regtype_ridge (l->lsys)) {	// Ridge penalty
				t[index] += alpha;
			} else {
				/* t = scale^2 * X(:,A)' * X(:,j) + scale^2 * lambda2 * J(:,A)' * J(:,j) */
				size_t			pj = linsys_get_pj (l->lsys);
				const double	*r = linsys_get_penalty (l->lsys);
				const double	*rj = r + LINSYS_INDEX_OF_MATRIX (0, j, pj);
				/* rtrj = J(:,A)' * J(:,j) */
				double			*rtrj = larsen_xa_transpose_dot_y (l, pj, 1., r, rj);
				/* t += scale^2 * lambda2 * J(:,A)' * J(:,j) */
				daxpy_ (LINSYS_CINTP (l->sizeA), &alpha, rtrj, &ione, t, &ione);
				free (rtrj);
			}
		}

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
	dcopy_ (LINSYS_CINTP (l->sizeA), s, &ione, l->w, &ione);

	/* cholesky update and solve equiangular equation
	 * TODO: create new methods to solve equiangular equation (using QR, SVD etc.)
	 * and switch them in the following line
	 */
	info = update_chol (l);
	if (!check_info ("update_chol", info)) return false;

	info = larsen_linalg_cholesky_svx (l->sizeA, l->chol, l->sizeA, l->w);
	if (!check_info ("cholesky_svx", info)) return false;
	/* end */

	l->absA = 1. / sqrt (ddot_ (LINSYS_CINTP (l->sizeA), s, &ione, l->w, &ione));
	free (s);

	/* w = (ZA' * ZA)^-1 * s(A) * absA */
	dscal_ (LINSYS_CINTP (l->sizeA), &l->absA, l->w, &ione);

	/* In exactlly
	 *
	 *     uA = ZA * w = scale * [XA * w; sqrt(lambda2) * EA * w],
	 *     dim(uA) = n + p (size of nonzero EA * w is |A|).
	 *
	 * But in this program, to save memory, only first n part of uA is stored
	 *
	 *     u = scale * XA * w,
	 *     dim(u) = n,
	 *
	 * because the later part of uA ( scale * sqrt(lambda2) * w ) is not appear
	 * in the positive on the LARS-EN algorithm.
	 */
	if (l->u) free (l->u);
	{
		size_t			n = linsys_get_n (l->lsys);
		double			scale = linsys_get_scale (l->lsys);
		const double	*x = linsys_get_x (l->lsys);
		l->u = larsen_xa_dot_ya (l, n, scale, x, l->w);
	}

	return true;
}

/* call equiangular vector updater
 * if another routine is implenented, swicth them in here */
bool
update_equiangular (larsen *l)
{
	return update_equiangular_larsen_cholesky (l);
}
