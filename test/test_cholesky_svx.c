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
test_cholesky_svx (void)
{
	size_t	size = 50;

	double *a;
	double *l;
	double *x;
	double *y;
	double *b;

	int		n = (int) size;
	int		nn = n * n;
	double	nrm;

	/* posdef symmetry matrix */
	a = posdef_symmetry_random_matrix (size);

	/* y = a * x */
	x = random_uniform_array (size, 1);
	y = (double *) malloc (size * sizeof (double));
	dgemv_ ("N", &n, &n, &done, a, &n, x, &ione, &dzero, y, &ione);
	free (x);

	/* cholesky decomposition */
	l = (double *) malloc (size * size * sizeof (double));
	dcopy_ (&nn, a, &ione, l, &ione);
	test_linalg_cholesky_decomp (size, l, size);

	/* x1 = a^-1 * y */
	x = (double *) malloc (size * sizeof (double));
	dcopy_ (&n, y, &ione, x, &ione);
	larsen_linalg_cholesky_svx (size, l, size, x);
	free (l);

	/* b = a * x1 */
	b = (double *) malloc (size * sizeof (double));
	dgemv_ ("N", &n, &n, &done, a, &n, x, &ione, &dzero, b, &ione);
	free (x);
	free (a);

	/* y = - b + y */
	daxpy_ (&n, &dmone, b, &ione, y, &ione);
	free (b);
	nrm = dnrm2_ (&n, y, &ione);
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
