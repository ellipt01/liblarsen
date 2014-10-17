/*
 * example.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <larsen.h>

void
output_solutionpath (int iter, larsen *l)
{
	int			i;
	char		fn[80];
	FILE		*fp;
	double		*beta = larsen_copy_beta (l, true);

	for (i = 0; i < l->beta->nz; i++) {

		sprintf (fn, "beta%03d.res", i);

		if (iter == 0) fp = fopen (fn, "w");
		else fp = fopen (fn, "aw");
		if (fp == NULL) continue;

		fprintf (fp, "%d\t%.4e\t%.4e\n", iter, l->lambda1, beta[i]);
		fclose (fp);
	}
	free (beta);
	return;
}

void
fprintf_beta (FILE *stream, int iter, larsen *l)
{
	int		i;
	double	*beta = larsen_copy_beta (l, true);

	fprintf (stream, "beta[%04d, t = %.2f] = \n", iter, l->lambda1);
	for (i = 0; i < l->beta->nz; i++) fprintf (stream, "%.4e\n", beta[i]);
	fprintf (stream, "]\n");
	free (beta);

	return;
}

