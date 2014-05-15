/*
 * data.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdlib.h>
#include <larsen.h>

#include "larsen_private.h"

/* Centering each column of matrix:
 * x(:, j) -> x(:, j) - mean(x(:, j)) */
double *
larsen_centering (const size_t size1, const size_t size2, double *x)
{
	int		i, j;
	double	*mean = (double *) malloc (size2 * sizeof (double));
	for (j = 0; j < size2; j++) {
		double	*xj = x + LARSEN_INDEX_OF_MATRIX (0, j, size1);
		mean[j] = 0.;
		for (i = 0; i < size1; i++) mean[j] += xj[i];
		mean[j] /= (double) size1;
		for (i = 0; i < size1; i++) xj[i] -= mean[j];
	}
	return mean;
}

/* Normalizing each column of matrix:
 * x(:, j) -> x(:, j) / norm(x(:, j)) */
double *
larsen_normalizing (const size_t size1, const size_t size2, double *x)
{
	int		j;
	double	*nrm = (double *) malloc (size2 * sizeof (double));
	for (j = 0; j < size2; j++) {
		double	alpha;
		double	*xj = x + LARSEN_INDEX_OF_MATRIX (0, j, size1);
		nrm[j] = dnrm2_ (CINTP (size1), xj, &ione);
		alpha = 1. / nrm[j];
		dscal_ (CINTP (size1), &alpha, xj, &ione);
	}
	return nrm;
}
