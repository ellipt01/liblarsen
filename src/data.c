/*
 * data.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <larsen.h>

/* centering vector: y -> y - mean(y) */
void
larsen_centering (const size_t size1, const size_t size2, double *x)
{
	int		i, j;
	for (j = 0; j < size2; j++) {
		double	*xj = x + INDEX_OF_MATRIX (0, j, size1);
		double	meanj = 0.;
		for (i = 0; i < size1; i++) meanj += xj[i];
		meanj /= (double) size1;
		for (i = 0; i < size1; i++) xj[i] -= meanj;
	}
	return;
}

/* normalizing each column of matrix:
 * x(:, j) -> x(:, j) / norm(x(:, j)) */
void
larsen_normalizing (const size_t size1, const size_t size2, double *x)
{
	int		j;
	for (j = 0; j < size2; j++) {
		double	*xj = x + INDEX_OF_MATRIX (0, j, size1);
		double	nrmj = cblas_dnrm2 (size1, xj, 1);
		cblas_dscal (size1, 1. / nrmj, xj, 1);
	}
	return;
}
