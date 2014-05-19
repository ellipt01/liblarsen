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

extern int		test_linalg_cholesky_decomp (const size_t size, double *a, const size_t lda);

bool
test_cholesky_delete (void)
{
	int		i;
	size_t	size = 50;
	size_t	index = 30;

	double *a;
	double *b;

	double *la;
	double *lb;

	int		n = (int) size;
	int		nn = n * n;
	int		mm = (int) (size - 1) * (size - 1);
	double	nrm;

	/* posdef symmetry matrix */
	b = posdef_symmetry_random_matrix (size);

	/* cholesky decomposition of matrix a */
	lb = (double *) malloc (size * size * sizeof (double));
	dcopy_ (&nn, b, &ione, lb, &ione);
	test_linalg_cholesky_decomp (size, lb, size);

	a = (double *) malloc ((size - 1) * (size - 1) * sizeof (double));
	{
		int		i, j, m, n;
		for (i = 0, m = 0; i < size; i++) {
			if (i == index) continue;
			for (j = 0, n = 0; j < size; j++) {
				if (j == index) continue;
				a[LINSYS_INDEX_OF_MATRIX (m, n, size - 1)] = b[LINSYS_INDEX_OF_MATRIX (i, j, size)];
				n++;
			}
			m++;
		}
	}
	free (b);

	/* cholesky decomposition */
	la = (double *) malloc ((size - 1) * (size - 1) * sizeof (double));
	dcopy_ (&mm, a, &ione, la, &ione);
	test_linalg_cholesky_decomp (size - 1, la, size - 1);
	for (i = 0; i < size - 1; i++) {
		int		j;
		for (j = 0; j < i; j++) la[LINSYS_INDEX_OF_MATRIX(i, j, size - 1)] = 0.;
	}
	free (a);

	larsen_linalg_cholesky_delete (size, &lb, index);
	for (i = 0; i < size - 1; i++) {
		int		j;
		for (j = 0; j < i; j++) lb[LINSYS_INDEX_OF_MATRIX(i, j, size - 1)] = 0.;
	}

	/* la = - lb + la */
	daxpy_ (&mm, &dmone, lb, &ione, la, &ione);
	free (lb);
	nrm = dnrm2_ (&mm, la, &ione);
	free (la);

	return (nrm < 1.e-8);
}

int
main (void)
{
	srand (time (NULL));
	fprintf (stdout, "%d\n", test_cholesky_delete () ? 1 : 0);
	return 0;
}
