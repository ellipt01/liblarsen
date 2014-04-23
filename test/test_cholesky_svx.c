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
test_cholesky_svx (void)
{
	size_t	size = 50;

	double *a;
	double *l;
	double *x;
	double *y;
	double *b;

	double	nrm;

	/* posdef symmetry matrix */
	a = posdef_symmetry_random_matrix (size);

	/* y = a * x */
	x = random_uniform_array (size, 1);
	y = (double *) malloc (size * sizeof (double));
	cblas_dgemv (CblasColMajor, CblasNoTrans, size, size, 1, a, size, x, 1, 0., y, 1);
	free (x);

	/* cholesky decomposition */
	l = (double *) malloc (size * size * sizeof (double));
	cblas_dcopy (size * size, a, 1, l, 1);
	larsen_linalg_cholesky_decomp (size, l, size);

	/* x1 = a^-1 * y */
	x = (double *) malloc (size * sizeof (double));
	cblas_dcopy (size, y, 1, x, 1);
	larsen_linalg_cholesky_svx (size, l, size, x);
	free (l);

	/* b = a * x1 */
	b = (double *) malloc (size * sizeof (double));
	cblas_dgemv (CblasColMajor, CblasNoTrans, size, size, 1, a, size, x, 1, 0., b, 1);
	free (x);
	free (a);

	/* y = - b + y */
	cblas_daxpy (size, -1., b, 1, y, 1);
	free (b);
	nrm = cblas_dnrm2 (size, y, 1);
	free (y);

	return (nrm < 1.e-8);
}

/* test_cholesky_svx */
int
main (void)
{
	srand (time (NULL));
	fprintf (stdout, "%d\n", test_cholesky_svx () ? 1 : 0);
	return 0;
}
