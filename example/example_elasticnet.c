/*
 * example_elasticnet.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <larsen.h>
#include "example.h"

void
example_elasticnet (cl_matrix *x, cl_vector *y, double start, double dt, double stop, double lambda2, double gamma, int maxiter)
{
	int			iter = 0;
	larsen		*l = larsen_alloc (start, lambda2, x, y);

	do {
		if (!larsen_elasticnet (l, maxiter)) break;
		output_solutionpath (iter++, l);
		fprintf (stdout, "%d : lambda1 = %f, bic(%.2f) = %f\n", iter, l->lambda1, gamma, larsen_eval_bic (l, gamma));
		larsen_increment_lambda1 (l, dt);
	} while (larsen_loop_continue (l, stop));

	larsen_free (l);

	return;
}
