/*
 * equiangular.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include "larsen.h"

/* return Xa = X(A) */
static cl_matrix *
extruct_xa (larsen *l)
{
	int			i;
	size_t		size1 = l->x->size1;
	size_t		size2 = l->A->size;
	cl_matrix	*xa;
	cl_vector	*col;

	if (size2 <= 0) return NULL;

	xa = cl_matrix_alloc (size1, size2);
	col = cl_vector_alloc (size1);
	for (i = 0; i < size2; i++) {
		int		j = cl_vector_int_get (l->A, i);
		cl_matrix_get_col (col, l->x, j);
		cl_matrix_set_col (xa, i, col);
	}
	cl_vector_free (col);

	return xa;
}

/* s_i = sign(c_i) */
static cl_vector *
update_sign (larsen *l)
{
	int			i;
	cl_vector	*s = cl_vector_alloc (l->A->size);
	for (i = 0; i < l->A->size; i++) {
		int		j = cl_vector_int_get (l->A, i);
		double	sign = (cl_vector_get (l->c, j) >= 0.) ? 1. : -1.;
		cl_vector_set (s, i, sign);
	}
	return s;
}

/* do cholinsert or choldelete of specified index */
static void
update_chol (larsen *l, cl_matrix *xa)
{
	int			index = l->oper.index;

	if (l->oper.action == ACTIVESET_ACTION_ADD) {
		/*** insert a predictor ***/
		int			column = l->oper.column;
		cl_vector	*t;
		cl_vector	*xi = cl_vector_alloc (l->x->size1);
		cl_matrix_get_col (xi, l->x, column);

		t = cl_matrix_transpose_dot_vector (xa, xi);
		if (l->do_scaling) {
			cl_vector_scale (t, pow (l->scale, 2.));
			cl_vector_set (t, index, cl_vector_get (t, index) + l->lambda2 * pow (l->scale, 2.));
		}
		cl_vector_free (xi);
		if (l->A->size == 1 && cl_matrix_is_empty (l->chol)) {
			l->chol = cl_matrix_alloc (1, 1);
			cl_matrix_set (l->chol, 0, 0, cl_vector_get (t, 0));
			cl_linalg_cholesky_decomp (l->chol);
		} else {
			cl_linalg_cholesky_insert (l->chol, index, t);
		}
		cl_vector_free (t);

	} else if (l->oper.action == ACTIVESET_ACTION_DROP) {
		/*** delete a predictor ***/
		cl_linalg_cholesky_delete (l->chol, index);
	}

	return;
}

/* update equiangular vector l->u, l->w and l->absA for current active set */
static bool
update_equiangular_larsen (larsen *l)
{
	cl_matrix	*xa;
	cl_vector	*s;

	if (l->A->size <= 0) return false;

	xa = extruct_xa (l);

	update_chol (l, xa);

	s = update_sign (l);

	if (l->w) cl_vector_free (l->w);
	l->w = cl_vector_alloc (s->size);
	cl_vector_memcpy (l->w, s);

	{
		int		info = cl_linalg_cholesky_svx (l->chol, l->w);
		if (info != 0) fprintf (stderr, "WARNING : cholesky_svx info = %d : matrix is not positive definite.\n", info);
	}

	l->absA = 1. / sqrt (cl_vector_dot_vector (s, l->w));
	cl_vector_free (s);
	cl_vector_scale (l->w, l->absA);

	if (l->u) cl_vector_free (l->u);
	l->u = cl_matrix_dot_vector (xa, l->w);
	if (l->do_scaling) cl_vector_scale (l->u, l->scale);
	cl_matrix_free (xa);

	return true;
}

/* another routine to calculate equiangular vector is implenented, swicth them */
bool
update_equiangular (larsen *l)
{
	return update_equiangular_larsen (l);
}
