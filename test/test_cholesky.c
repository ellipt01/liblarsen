/*
 * test_cholesky.c
 *
 *  Created on: 2014/04/23
 *      Author: utsugi
 */

#include <stdlib.h>
#include <larsen.h>

#include "test_larsen.h"

const double	done = 1.;

#ifndef HAVE_LAPACK_H
extern void	dpotrf_ (char *uplo, const int *n, double *a, int *lda, int *info);
#endif

int
test_linalg_cholesky_decomp (const size_t size, double *a, const size_t lda)
{
	int		info;
	char	uplo;
	int		n;
	int		_lda;

	if (!a) linreg_error ("linreg_linalg_cholesky_decomp", "matrix is empty.", __FILE__, __LINE__);

	uplo = 'U';
	n = (int) size;
	_lda = (int) lda;
	dpotrf_ (&uplo, &n, a, &_lda, &info);
	return info;
}

double *
random_uniform_array (const size_t size1, const size_t size2)
{
	int		i;
	size_t	n = size1 * size2;
	double	*array = (double *) malloc (n * sizeof (double));

	for (i = 0; i < n; i++) array[i] = (double) rand () / RAND_MAX;

	return array;
}

double *
posdef_symmetry_random_matrix (const size_t size)
{
	int		i;
	int		n = (int) size;
	double	*x = (double *) malloc (size * size * sizeof (double));
	double	*x0 = random_uniform_array (size, size);
	dgemm_ ("T", "N", &n, &n, &n, &done, x0, &n, x0, &n, &dzero, x, &n);
	free (x0);
	for (i = 0; i < size; i++) x[i + i * size] += 1.;
	return x;
}
