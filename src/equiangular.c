/*
 * equiangular.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <larsen.h>

/* return XA = X(A) */
static c_matrix *
extruct_xa (larsen *l)
{
	int			i;
	size_t		size1 = l->x->size1;
	size_t		size2 = l->sizeA;
	c_vector	*xj;
	c_matrix	*xa;

	if (size2 <= 0) return NULL;

	xj = c_vector_alloc (size1);
	xa = c_matrix_alloc (size1, size2);
	for (i = 0; i < size2; i++) {
		int		j = l->A[i];
		c_matrix_get_col (xj, l->x, j);
		c_matrix_set_col (xa, i, xj);
	}
	c_vector_free (xj);

	return xa;
}

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

/* do cholinsert or choldelete of specified index */
static int
update_chol (larsen *l, c_matrix *xa)
{
	int		info = 0;
	int		index = l->oper.index;

	if (l->oper.action == ACTIVESET_ACTION_ADD) {
		/*** insert a predictor ***/
		int			column = l->oper.column;
		c_vector	*t;
		c_vector	*xi = c_vector_alloc (l->x->size1);
		c_matrix_get_col (xi, l->x, column);
		t = c_matrix_transpose_dot_vector (l->scale2, xa, xi, 0.);
		c_vector_free (xi);
		if (l->is_elnet) c_vector_set (t, index, c_vector_get (t, index) + l->lambda2 * l->scale2);
		if (l->sizeA == 1 && c_matrix_is_empty (l->chol)) {
			l->chol = c_matrix_alloc (1, 1);
			c_matrix_set (l->chol, 0, 0, c_vector_get (t, 0));
			info = c_linalg_cholesky_decomp (1, l->chol->data);
		} else {
			fprintf (stderr, "l->sizeA = %d\n", (int) l->sizeA);
			fprintf (stderr, "0: l->chol[%d, %d]\n", (int) l->chol->size1, (int) l->chol->size2);
			c_matrix_fprintf (stderr, l->chol, "%.3f");
			info = c_linalg_cholesky_insert (l->sizeA - 1, l->chol->data, index, t->data);
			fprintf (stderr, "1: l->chol[%d, %d]\n", (int) l->chol->size1, (int) l->chol->size2);
			l->chol->size1++;
			l->chol->size2++;
			fprintf (stderr, "2: l->chol[%d, %d]\n", (int) l->chol->size1, (int) l->chol->size2);
			c_matrix_fprintf (stderr, l->chol, "%.3f");
	}
		c_vector_free (t);

	} else if (l->oper.action == ACTIVESET_ACTION_DROP) {
		/*** delete a predictor ***/
		c_linalg_cholesky_delete (l->sizeA + 1, l->chol->data, index);
	}

	return info;
}

/* update equiangular vector and its relevant (l->u, l->w and l->absA)
 * for current active set using cholinsert / delete routine */
static bool
update_equiangular_larsen_cholesky (larsen *l)
{
	int			info;
	c_matrix	*xa;
	double		*s;

	if (l->sizeA <= 0) return false;

	xa = extruct_xa (l);
	
	s = update_sign (l);

	if (l->w) free (l->w);
	l->w = (double *) malloc (l->sizeA * sizeof (double));
	cblas_dcopy (l->sizeA, s, 1, l->w, 1);

	/* cholesky update and solve equiangular equation
	 * TODO: create new method and move the following lines to it
	 * to be able to switch some methods
	 * */
	info = update_chol (l, xa);
	if (!check_info ("update_chol", info)) return false;

	info = c_linalg_cholesky_svx (l->sizeA, l->chol->data, l->w);
	if (!check_info ("cholesky_svx", info)) return false;
	/* end */

	l->absA = 1. / sqrt (cblas_ddot (l->sizeA, s, 1, l->w, 1));
	free (s);
	cblas_dscal (l->sizeA, l->absA, l->w, 1);

	if (!l->u) l->u = (double *) malloc (l->n * sizeof (double));
	cblas_dgemv (CblasColMajor, CblasNoTrans, l->n, l->sizeA, l->scale, xa->data, l->n, l->w, 1, 0., l->u, 1);

	c_matrix_free (xa);

	return true;
}

/* call equiangular vector updater
 * if another routine is implenented, swicth them in here */
bool
update_equiangular (larsen *l)
{
	return update_equiangular_larsen_cholesky (l);
}
