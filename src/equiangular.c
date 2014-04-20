/*
 * equiangular.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdlib.h>
#include <math.h>
#include <larsen.h>

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
	bool	is_valid = false;

	if (info == 0) is_valid = true;
	else if (info < 0) fprintf (stderr, "WARNING: %s : the %d-th argument had an illegal value.\n", func_name, - info - 1);
	else if (info > 0) {
		fprintf (stderr, "WARNING: %s : the leading minor of order %d is not positive definite.\n", func_name, info);
		fprintf (stderr, "         factorization could not be completed.\n");
	}
	return is_valid;
}

/* y = alpha * XA * z */
static double *
xa_dot_y (larsen *l, double alpha, double *z)
{
	int		j;
	double	*y = (double *) malloc (l->n * sizeof (double));
	double	*zp = (double *) malloc (l->p * sizeof (double));
	/* zp(j) = z(j), j in A, els = 0, j not in A */
	for (j = 0; j < l->p; j++) zp[j] = 0.;
	for (j = 0; j < l->sizeA; j++) zp[l->A[j]] = z[j];
	/* y = X * zp = X(A, :) * z */
	cblas_dgemv (CblasColMajor, CblasNoTrans, l->n, l->p, alpha, l->x, l->n, zp, 1, 0., y, 1);
	free (zp);

	return y;
}

/* y = alpha * XA' * z */
static double *
xa_transpose_dot_y (larsen *l, const double alpha, const double *z)
{
	int		j;
	double	*y = (double *) malloc (l->sizeA * sizeof (double));

	/* another version with cblas_dgemv */
	/*	double	*yp = (double *) malloc (l->p * sizeof (double));
		cblas_dgemv (CblasColMajor, CblasTrans, l->n, l->p, alpha, l->x, l->n, z, 1, 0., yp, 1);
		for (j = 0; j < l->sizeA; j++) y[j] = yp[l->A[j]];
		free (yp); */

	/* The following is more fast when l->sizeA is not huge */
	for (j = 0; j < l->sizeA; j++) {
		const double	*xaj = l->x + INDEX_OF_MATRIX (0, l->A[j], l->n);
		/* y[j] = alpha * X(:, A[j])' * z */
		y[j] = alpha * cblas_ddot (l->n, xaj, 1, z, 1);
	}
	return y;
}

/* do cholinsert or choldelete of specified index */
static int
update_chol (larsen *l)
{
	int		info = 0;
	int		index = l->oper.index_of_A;

	if (l->oper.action == ACTIVESET_ACTION_ADD) {
		/*** insert a predictor ***/
		int				j = l->oper.column_of_X;
		double			*t = (double *) malloc (l->sizeA * sizeof (double));
		const double	*xj = l->x + INDEX_OF_MATRIX (0, j, l->n);

		/* t = scale^2 * X(:,A)' * X(:,j) */
		t = xa_transpose_dot_y (l, l->scale2, xj);
		/* t += lambda2 * scale^2 */
		if (l->is_elnet) t[index] += l->lambda2 * l->scale2;
		info = larsen_linalg_cholesky_insert (l->sizeA - 1, &l->chol, index, t);
		free (t);

	} else if (l->oper.action == ACTIVESET_ACTION_DROP) {
		/*** delete a predictor ***/
		larsen_linalg_cholesky_delete (l->sizeA + 1, &l->chol, index);
	}

	return info;
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

	info = larsen_linalg_cholesky_svx (l->sizeA, l->chol, l->sizeA, l->w);
	if (!check_info ("cholesky_svx", info)) return false;
	/* end */

	l->absA = 1. / sqrt (cblas_ddot (l->sizeA, s, 1, l->w, 1));
	free (s);
	cblas_dscal (l->sizeA, l->absA, l->w, 1);

	/* u = scale * XA * w */
	if (l->u) free (l->u);
	l->u = xa_dot_y (l, l->scale, l->w);

	return true;
}

/* call equiangular vector updater
 * if another routine is implenented, swicth them in here */
bool
update_equiangular (larsen *l)
{
	return update_equiangular_larsen_cholesky (l);
}
