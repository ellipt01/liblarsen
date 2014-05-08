/*
 * test_cholesky.c
 *
 *  Created on: 2014/04/23
 *      Author: utsugi
 */

#include <stdlib.h>
#include <larsen.h>

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

	if (!a) larsen_linalg_error ("larsen_linalg_cholesky_decomp", "matrix is empty.");

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
	double	*x = (double *) malloc (size * size * sizeof (double));
	double	*x0 = random_uniform_array (size, size);
	cblas_dgemm (CblasColMajor, CblasTrans, CblasNoTrans, size, size, size, 1., x0, size, x0, size, 0., x, size);
	free (x0);
	for (i = 0; i < size; i++) x[i + i * size] += 1.;
	return x;
}