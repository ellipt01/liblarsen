/*
 * utils.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <larsen.h>

larsen *
larsen_alloc (double lambda1, double lambda2, const c_matrix *x, const c_vector *y)
{
	int		i;
	larsen	*l = (larsen *) malloc (sizeof (larsen));

	l->stop_loop = false;

	l->lambda1 = lambda1;
	l->lambda2 = lambda2;

	/* vector and matrix view of original y and x */
	l->x = c_matrix_view_array (x->size1, x->size2, x->lda, x->data);
	l->y = c_vector_view_array (y->size, y->stride, y->data);

	l->scale = 1.;
	l->do_scaling = (lambda2 > DBL_EPSILON);
	if (l->do_scaling) l->scale /= sqrt (1. + lambda2);

	/* correlation */
	l->sup_c = 0.;
	l->c = NULL;

	/* active set */
	l->oper.action = ACTIVESET_ACTION_ADD;
	l->oper.column = -1;
	l->oper.index = -1;

	l->A = c_vector_int_alloc (x->size2);
	l->A->size = 0;
	l->Ac = c_vector_int_alloc (x->size2);
	for (i = 0; i < x->size2; i++) c_vector_int_set (l->Ac, i, i);

	/* active set */
	l->absA = 0.;
	l->u = NULL;
	l->w = NULL;

	/* solution */
	l->beta = c_vector_alloc (x->size2);
	c_vector_set_zero (l->beta);
	l->mu = c_vector_alloc (y->size);
	c_vector_set_zero (l->mu);

	/* backup of solution */
	l->beta_prev = c_vector_alloc (x->size2);
	c_vector_set_zero (l->beta_prev);
	l->mu_prev = c_vector_alloc (y->size);
	c_vector_set_zero (l->mu_prev);

	/* interpolation */
	l->interp = false;
	l->stepsize_intr = 0.;
	l->beta_intr = c_vector_alloc (x->size2);
	l->mu_intr = c_vector_alloc (y->size);

	/* cholesky factorization */
	l->chol = NULL;

	return l;
}

void
larsen_free (larsen *l)
{
	if (l) {
		if (!c_vector_is_empty (l->c)) c_vector_free (l->c);


		if (!c_vector_int_is_empty (l->A)) c_vector_int_free (l->A);
		if (!c_vector_int_is_empty (l->Ac)) c_vector_int_free (l->Ac);

		if (!c_matrix_is_empty (l->x)) c_matrix_free ((c_matrix *) l->x);
		if (!c_vector_is_empty (l->y)) c_vector_free ((c_vector *) l->y);

		if (!c_vector_is_empty (l->u)) c_vector_free (l->u);
		if (!c_vector_is_empty (l->w)) c_vector_free (l->w);

		if (!c_vector_is_empty (l->beta)) c_vector_free (l->beta);
		if (!c_vector_is_empty (l->mu)) c_vector_free (l->mu);

		if (!c_vector_is_empty (l->beta_prev)) c_vector_free (l->beta_prev);
		if (!c_vector_is_empty (l->mu_prev)) c_vector_free (l->mu_prev);

		if (!c_vector_is_empty (l->beta_intr)) c_vector_free (l->beta_intr);
		if (!c_vector_is_empty (l->mu_intr)) c_vector_free (l->mu_intr);

		if (!c_matrix_is_empty (l->chol)) c_matrix_free (l->chol);

		free (l);
	}
	return;
}
