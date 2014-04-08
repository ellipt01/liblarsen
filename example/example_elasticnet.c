/*
 * example_elasticnet.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <larsen.h>
#include "example.h"

/*
 * Calculate lasso or elastic net solution path for lambda1 in [start : dt : stop].
 * The value of extended BIC (Chen and Chen, 2008) are also calculated for each lambda1.
 */
void
example_elasticnet (c_matrix *x, c_vector *y, double start, double dt, double stop, double lambda2, double gamma, int maxiter)
{
	int			iter = 0;
	double		t = start;
	larsen		*l = larsen_alloc (t, lambda2, x, y);

	while (larsen_elasticnet (l, maxiter)) {
		output_solutionpath (iter++, l);
		fprintf (stdout, "%d : lambda1 = %f, bic(%.2f) = %f\n", iter, l->lambda1, gamma, larsen_eval_bic (l, gamma));
		t += dt;
		if (t > stop) break;
		larsen_set_lambda1 (l, t);
	}

	larsen_free (l);

	return;
}
