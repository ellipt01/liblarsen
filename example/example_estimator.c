/*
 * example_elasticnet.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdio.h>
#include <larsen.h>
#include "example.h"

/*
 * Calculate solution path of L-1 or L-1 and L-n competitive penalized regression [start : dt : stop].
 * The value of extended BIC (Chen and Chen, 2008) are also calculated for each lambda1.
 */
void
example_estimator (const linreg *lreg, double start, double dt, double stop, double gamma, int maxiter)
{
	int			iter = 0;
	double		t = start;
	larsen		*l = larsen_alloc (lreg, t);

	if (l == NULL) return;

	while (larsen_estimator (l, maxiter)) {
		output_solutionpath (iter++, l);
		fprintf (stdout, "%d : lambda1 = %f, bic(%.2f) = %f", iter, larsen_get_lambda1 (l, false), gamma, larsen_eval_bic (l, gamma));
		fprintf (stdout, ", df = %d\n", l->sizeA);
		t += dt;
		if (t > stop) break;
		larsen_set_lambda1 (l, t);
	}

	larsen_free (l);

	return;
}
