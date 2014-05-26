/*
 * example_cdescent.c
 *
 *  Created on: 2014/05/27
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <cdescent.h>

static void
output_solutionpath_cdescent (int iter, const cdescent *cd)
{
	int			i;
	char		fn[80];
	FILE		*fp;
	size_t		p = cd->lreg->p;
	double		nrm1 = cd->nrm1 / cd->lreg->scale2;
	double		*beta = cdescent_copy_beta (cd, true);

	for (i = 0; i < p; i++) {

		sprintf (fn, "beta%03d.res", i);

		if (iter == 0) fp = fopen (fn, "w");
		else fp = fopen (fn, "aw");
		if (fp == NULL) continue;

		fprintf (fp, "%d\t%.4e\t%.4e\n", iter, nrm1, beta[i]);
		fclose (fp);
	}
	free (beta);
	return;
}

void
example_cdescent_cyclick (const linreg *lreg, double start, double dt, double stop, double tol, int maxiter)
{
	int			iter = 0;
	double		t = start;
	cdescent	*cd = cdescent_alloc (lreg, t, tol);

	while (t <= stop) {

		if (!cdescent_cyclick (cd, maxiter)) break;
		output_solutionpath_cdescent (iter++, cd);

		t += dt;
		cdescent_set_lambda1 (cd, t);
	}

	return;
}
