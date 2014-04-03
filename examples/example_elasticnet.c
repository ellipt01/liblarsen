/*
 * example_elasticnet.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include "larsen.h"
#include "examples/examples.h"

void
dummy_ (void)
{

	return;
}

void
example_elasticnet (cl_matrix *x, cl_vector *y, double start, double dt, double stop, double lambda2, int maxiter)
{
	int			iter = 0;
	larsen		*l = larsen_alloc (start, lambda2, x, y);

	do {
		if (!larsen_elasticnet (l, maxiter)) break;
		output_solutionpath (iter++, l);
//		fprintf (stdout, "%d : lambda1 = %f, bic = %f\n", iter, l->lambda1, larsen_eval_bic (l, 0.));
		fprintf (stdout, "%d : lambda1 = %f, bic = %f\n", iter, l->lambda1, larsen_eval_bic (l, 0.5));
		larsen_increment_lambda1 (l, dt);
		if (l->lambda1 > 2861.) dummy_ ();
	} while (larsen_loop_continue (l, stop));

	larsen_free (l);

	return;
}
