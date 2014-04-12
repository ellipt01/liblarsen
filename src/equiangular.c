/*
 * equiangular.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdlib.h>
#include <math.h>
#include <larsen.h>

/* s_i = sign(c_i) */
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

/* check result of cholesky decomp, insert or delete */
static bool
check_info (const char *func_name, const int info)
{
	bool	is_valid = false;

	if (info == 0) is_valid = true;
	else if (info < 0) fprintf (stderr, "WARNING: %s : the %d-th argument had an illegal value.\n", func_name, - info - 1);
	else if (info > 0) {
		fprintf (stderr, "WARNING: %s : the leading minor of order %d is not positive definite.\n", func_name, info);
		fprintf (stderr, "         factorization could not be completed.\n");
	}
	return is_valid;
}

/* t = alpha * XA' * x */
static double *
xa_transpose_dot_xj (larsen *l, double alpha, const double *xj)
{
	int		j;
	double	*t = (double *) malloc (l->sizeA * sizeof (double));
	for (j = 0; j < l->sizeA; j++) {
		int				k = l->A[j];
		const double	*xk = l->x + index_of_matrix (0, k, l->n);
		/* y[j] = alpha * X(:, A[j])' * xj */
		t[j] = alpha * cblas_ddot (l->n, xk, 1, xj, 1);
	}
	return t;
}

/* do cholinsert or choldelete of specified index */
static int
update_chol (larsen *l)
{
	int		info = 0;
	int		index = l->oper.index;

	if (l->oper.action == ACTIVESET_ACTION_ADD) {
		/*** insert a predictor ***/
		int				j = l->oper.column;
		double			*t = (double *) malloc (l->sizeA * sizeof (double));
		const double	*xj = l->x + index_of_matrix (0, j, l->n);

		/* t = scale^2 * X(:,A)' * X(:,j) */
		t = xa_transpose_dot_xj (l, l->scale2, xj);
		/* t += lambda2 * scale^2 */
		if (l->is_elnet) t[index] += l->lambda2 * l->scale2;

		if (l->sizeA == 1 && l->chol == NULL) {
			l->chol = (double *) malloc (1 * sizeof (double));
			l->chol[0] = t[0];
			info = clinalg_cholesky_decomp (1, l->chol, 1);
		} else {
			info = clinalg_cholesky_insert (l->sizeA - 1, &l->chol, index, t);
		}
		free (t);

	} else if (l->oper.action == ACTIVESET_ACTION_DROP) {
		/*** delete a predictor ***/
		clinalg_cholesky_delete (l->sizeA + 1, &l->chol, index);
	}

	return info;
}

/* u = scale * XA * w */
static double *
xa_dot_w (larsen *l, double alpha, double *w)
{
	int		i, j;
	double	*u = (double *) malloc (l->n * sizeof (double));
	double	*row = (double *) malloc (l->sizeA * sizeof (double));
	for (i = 0; i < l->n; i++) {
		for (j = 0; j < l->sizeA; j++) {
			int		k = l->A[j];
			row[j] = l->x[index_of_matrix (i, k, l->n)];
		}
		/* u[i] = alpha * X(:, A) * w */
		u[i] = alpha * cblas_ddot (l->sizeA, row, 1, l->w, 1);
	}
	free (row);
	return u;
}

/* update equiangular vector and its relevant (l->u, l->w and l->absA)
 * for current active set using cholinsert / delete routine */
static bool
update_equiangular_larsen_cholesky (larsen *l)
{
	int			info;
	double		*s;

	if (l->sizeA <= 0) return false;

	s = update_sign (l);

	if (l->w) free (l->w);
	l->w = (double *) malloc (l->sizeA * sizeof (double));
	cblas_dcopy (l->sizeA, s, 1, l->w, 1);

	/* cholesky update and solve equiangular equation
	 * TODO: create new method and move the following lines to it
	 * to be able to switch some methods
	 * */
	info = update_chol (l);
	if (!check_info ("update_chol", info)) return false;

	info = clinalg_cholesky_svx (l->sizeA, l->chol, l->sizeA, l->w);
	if (!check_info ("cholesky_svx", info)) return false;
	/* end */

	l->absA = 1. / sqrt (cblas_ddot (l->sizeA, s, 1, l->w, 1));
	free (s);
	cblas_dscal (l->sizeA, l->absA, l->w, 1);

	/* u = scale * XA * w */
	if (l->u) free (l->u);
	l->u = xa_dot_w (l, l->scale, l->w);

	return true;
}

/* call equiangular vector updater
 * if another routine is implenented, swicth them in here */
bool
update_equiangular (larsen *l)
{
	return update_equiangular_larsen_cholesky (l);
}
