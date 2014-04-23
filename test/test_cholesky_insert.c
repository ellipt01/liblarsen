/*
 * test_cholesky.c
 *
 *  Created on: 2014/04/23
 *      Author: utsugi
 */

#include <stdlib.h>
#include <stdbool.h>
#include <larsen.h>

#include "test_larsen.h"

bool
test_cholesky_insert (void)
{
	int		i;
	size_t	size = 50;
	size_t	index = size - 1;

	double *a;
	double *b;
	double *u;

	double *la;
	double *lb;

	double	nrm;

	/* posdef symmetry matrix */
	a = posdef_symmetry_random_matrix (size);

	/* cholesky decomposition of matrix a */
	la = (double *) malloc (size * size * sizeof (double));
	cblas_dcopy (size * size, a, 1, la, 1);
	larsen_linalg_cholesky_decomp (size, la, size);
	for (i = 0; i < size; i++) {
		int		j;
		for (j = 0; j < i; j++) la[INDEX_OF_MATRIX(i, j, size)] = 0.;
	}

	/* b = a(1:end-1, 1:end-1), u = a(:, index) */
	b = (double *) malloc ((size - 1) * (size - 1) * sizeof (double));
	u = (double *) malloc (size * sizeof (double));
	{
		int		i, j;
		for (i = 0; i < size - 1; i++) {
			for (j = 0; j < size - 1; j++) {
				b[INDEX_OF_MATRIX (i, j, size - 1)] = a[INDEX_OF_MATRIX (i, j, size)];
			}
		}
		for (i = 0; i < size; i++) u[i] = a[INDEX_OF_MATRIX (i, index, size)];
	}
	free (a);

	/* cholesky decomposition */
	lb = (double *) malloc ((size - 1) * (size - 1) * sizeof (double));
	cblas_dcopy ((size - 1) * (size - 1), b, 1, lb, 1);
	larsen_linalg_cholesky_decomp (size - 1, lb, size - 1);
	free (b);

	larsen_linalg_cholesky_insert (size - 1, &lb, index, u);
	free (u);
	for (i = 0; i < size; i++) {
		int		j;
		for (j = 0; j < i; j++) lb[INDEX_OF_MATRIX(i, j, size)] = 0.;
	}

	/* la = - lb + la */
	cblas_daxpy (size * size, -1., lb, 1, la, 1);
	free (lb);
	nrm = cblas_dnrm2 (size * size, la, 1);
	free (la);

	return (nrm < 1.e-8);
}

int
main (void)
{
	srand (time (NULL));
	fprintf (stdout, "%d\n", test_cholesky_insert () ? 1 : 0);
	return 0;
}
