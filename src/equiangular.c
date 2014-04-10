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
	size_t		size2 = l->A->size;
	c_vector	*xj;
	c_matrix	*xa;

	if (size2 <= 0) return NULL;

	xj = c_vector_alloc (size1);
	xa = c_matrix_alloc (size1, size2);
	for (i = 0; i < size2; i++) {
		int		j = c_vector_int_get (l->A, i);
		c_matrix_get_col (xj, l->x, j);
		c_matrix_set_col (xa, i, xj);
	}
	c_vector_free (xj);

	return xa;
}

/* s_i = sign(c_i) */
static c_vector *
update_sign (larsen *l)
{
	int			i;
	c_vector	*s = c_vector_alloc (l->A->size);
	for (i = 0; i < l->A->size; i++) {
		int		j = c_vector_int_get (l->A, i);
		double	sign = (c_vector_get (l->c, j) >= 0.) ? 1. : -1.;
		c_vector_set (s, i, sign);
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
		if (l->A->size == 1 && c_matrix_is_empty (l->chol)) {
			l->chol = c_matrix_alloc (1, 1);
			c_matrix_set (l->chol, 0, 0, c_vector_get (t, 0));
			info = c_linalg_cholesky_decomp (l->chol);
		} else {
			info = c_linalg_cholesky_insert (l->chol, index, t);
		}
		c_vector_free (t);

	} else if (l->oper.action == ACTIVESET_ACTION_DROP) {
		/*** delete a predictor ***/
		c_linalg_cholesky_delete (l->chol, index);
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
	c_vector	*s;

	if (l->A->size <= 0) return false;

	xa = extruct_xa (l);
	
	s = update_sign (l);

	if (!c_vector_is_empty (l->w)) c_vector_free (l->w);
	l->w = c_vector_alloc (s->size);
	c_vector_memcpy (l->w, s);

	/* cholesky update and solve equiangular equation
	 * TODO: create new method and move the following lines to it
	 * to be able to switch some methods
	 * */
	info = update_chol (l, xa);
	if (!check_info ("update_chol", info)) return false;

	info = c_linalg_cholesky_svx (l->chol, l->w);
	if (!check_info ("cholesky_svx", info)) return false;
	/* end */

	l->absA = 1. / sqrt (c_vector_dot_vector (s, l->w));
	c_vector_free (s);
	c_vector_scale (l->w, l->absA);

	if (!c_vector_is_empty (l->u)) c_vector_free (l->u);

	l->u = c_matrix_dot_vector (l->scale, xa, l->w, 0.);

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
