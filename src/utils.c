/*
 * utils.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include "larsen.h"

larsen *
larsen_alloc (double lambda1, double lambda2, const cl_matrix *x, const cl_vector *y)
{
	int		i;
	larsen	*l = (larsen *) malloc (sizeof (larsen));

	l->stop_loop = false;

	l->lambda1 = lambda1;
	l->lambda2 = lambda2;
	l->x = x;
	l->y = y;

	l->scale = 1.;
	l->do_scaling = (lambda2 > CL_DBL_EPSILON);
	if (l->do_scaling) l->scale /= sqrt (1. + lambda2);

	/* correlation */
	l->sup_c = 0.;
	l->c = NULL;

	/* active set */
	l->oper.action = ACTIVESET_ACTION_ADD;
	l->oper.column = -1;
	l->oper.index = -1;

	l->A = cl_vector_int_alloc (x->size2);
	l->A->size = 0;
	l->Ac = cl_vector_int_alloc (x->size2);
	for (i = 0; i < x->size2; i++) cl_vector_int_set (l->Ac, i, i);

	/* active set */
	l->absA = 0.;
	l->u = NULL;
	l->w = NULL;

	/* solution */
	l->beta = cl_vector_alloc (x->size2);
	cl_vector_set_zero (l->beta);
	l->mu = cl_vector_alloc (y->size);
	cl_vector_set_zero (l->mu);

	/* backup of solution */
	l->beta_prev = cl_vector_alloc (x->size2);
	cl_vector_set_zero (l->beta_prev);
	l->mu_prev = cl_vector_alloc (y->size);
	cl_vector_set_zero (l->mu_prev);

	/* interpolation */
	l->interp = false;
	l->stepsize_intr = 0.;
	l->beta_intr = cl_vector_alloc (x->size2);
	l->mu_intr = cl_vector_alloc (y->size);

	/* cholesky factorization */
	l->chol = NULL;

	return l;
}

void
larsen_free (larsen *l)
{
	if (l) {
		if (!cl_vector_is_empty (l->c)) cl_vector_free (l->c);

		if (!cl_vector_int_is_empty (l->A)) cl_vector_int_free (l->A);
		if (!cl_vector_int_is_empty (l->Ac)) cl_vector_int_free (l->Ac);

		if (!cl_vector_is_empty (l->u)) cl_vector_free (l->u);
		if (!cl_vector_is_empty (l->w)) cl_vector_free (l->w);

		if (!cl_vector_is_empty (l->beta)) cl_vector_free (l->beta);
		if (!cl_vector_is_empty (l->mu)) cl_vector_free (l->mu);

		if (!cl_vector_is_empty (l->beta_prev)) cl_vector_free (l->beta_prev);
		if (!cl_vector_is_empty (l->mu_prev)) cl_vector_free (l->mu_prev);

		if (!cl_vector_is_empty (l->beta_intr)) cl_vector_free (l->beta_intr);
		if (!cl_vector_is_empty (l->mu_intr)) cl_vector_free (l->mu_intr);

		if (!cl_matrix_is_empty (l->chol)) cl_matrix_free (l->chol);

		free (l);
	}
	return;
}
