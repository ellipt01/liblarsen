/*
 * data.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdlib.h>
#include <larsen.h>

/* centering vector: y -> y - mean(y) */
double *
larsen_centering (const size_t size1, const size_t size2, double *x)
{
	int		i, j;
	double	*mean = (double *) malloc (size2 * sizeof (double));
	for (j = 0; j < size2; j++) {
		double	*xj = x + INDEX_OF_MATRIX (0, j, size1);
		mean[j] = 0.;
		for (i = 0; i < size1; i++) mean[j] += xj[i];
		mean[j] /= (double) size1;
		for (i = 0; i < size1; i++) xj[i] -= mean[j];
	}
	return mean;
}

/* normalizing each column of matrix:
 * x(:, j) -> x(:, j) / norm(x(:, j)) */
double *
larsen_normalizing (const size_t size1, const size_t size2, double *x)
{
	int		j;
	double	*nrm = (double *) malloc (size2 * sizeof (double));
	for (j = 0; j < size2; j++) {
		double	*xj = x + INDEX_OF_MATRIX (0, j, size1);
		nrm[j] = cblas_dnrm2 (size1, xj, 1);
		cblas_dscal (size1, 1. / nrm[j], xj, 1);
	}
	return nrm;
}
