/*
 * example_elasticnet.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <larsen.h>
#include "example.h"

static void
output_solutionpath (int iter, larsen *l)
{
	int			i;
	char		fn[80];
	FILE		*fp;
	size_t		p = l->lreg->p;
	double		lambda1 = larsen_get_lambda1 (l, false);
	double		*beta = larsen_copy_beta (l, true);

	for (i = 0; i < p; i++) {

		sprintf (fn, "beta%03d.res", i);

		if (iter == 0) fp = fopen (fn, "w");
		else fp = fopen (fn, "aw");
		if (fp == NULL) continue;

		fprintf (fp, "%d\t%.4e\t%.4e\n", iter, lambda1, beta[i]);
		fclose (fp);
	}
	free (beta);
	return;
}

/*
 * Calculate solution path of L-1 or L-1 and L-n competitive penalized regression [start : dt : stop].
 * The value of extended BIC (Chen and Chen, 2008) are also calculated for each lambda1.
 */
void
example_larsen (const linreg *lreg, double start, double dt, double stop, double gamma, int maxiter)
{
	int			iter = 0;
	double		t = start;
	larsen		*l = larsen_alloc (lreg, t);

	if (l == NULL) return;

	while (larsen_estimater (l, maxiter)) {
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
