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
test_cholesky_decomp (void)
{
	int		i;
	size_t	size = 50;

	double *a;
	double	*l;
	double	*ltl;

	double	nrm;

	/* posdef symmetry matrix */
	a = posdef_symmetry_random_matrix (size);

	/* cholesky decomposition */
	l = (double *) malloc (size * size * sizeof (double));
	cblas_dcopy (size * size, a, 1, l, 1);
	larsen_linalg_cholesky_decomp (size, l, size);
	for (i = 0; i < size; i++) {
		int		j;
		for (j = 0; j < i; j++) l[INDEX_OF_MATRIX(i, j, size)] = 0.;
	}

	/* l^T * l */
	ltl = (double *) malloc (size * size * sizeof (double));
	cblas_dgemm (CblasColMajor, CblasTrans, CblasNoTrans, size, size, size, 1., l, size, l, size, 0., ltl, size);
	free (l);

	/* a - ltl */
	cblas_daxpy (size * size, -1., ltl, 1, a, 1);
	free (ltl);

	nrm = cblas_dnrm2 (size * size, a, 1);
	free (a);

	return (nrm < 1.e-8);
}

/* test_cholesky_decomp */
int
main (void)
{
	srand (time (NULL));
	fprintf (stdout, "%d\n", test_cholesky_decomp () ? 1 : 0);
	return 0;
}
