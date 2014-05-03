/*
 * larsen_linalg_cholesky.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

//#include <stdio.h>
#include <stdlib.h>
#include <larsen_linalg.h>

extern void	larsen_linalg_error (const char * function_name, const char *error_msg);
extern void	matrix_add_rowcols (const size_t m, const size_t n, double **a, const size_t dm, const size_t dn);
extern void	matrix_remove_rowcols (const size_t m, const size_t n, double **a, const size_t dm, const size_t dn);

int
larsen_linalg_cholesky_decomp (const size_t size, double *a, const size_t lda)
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

int
larsen_linalg_cholesky_svx (const size_t size, double *l, const size_t lda, double *b)
{
	int		info;
	char	uplo;
	int		n;
	int		nrhs;
	int		_lda;

	if (!l) larsen_linalg_error ("larsen_linalg_cholesky_svx", "first matrix is empty.");
	if (!b) larsen_linalg_error ("larsen_linalg_cholesky_svx", "second matrix is empty.");

	uplo = 'U';
	n = (int) size;
	nrhs = 1;
	_lda = (int) lda;
	dpotrs_ (&uplo, &n, &nrhs, l, &_lda, b, &n, &info);

	return info;
}

int
larsen_linalg_cholesky_insert (const size_t size, double **r, const int index, double *u)
{
	int			info;
	int			n;
	int			ldr;
	int			j;
	double		*w;

	if (index < 0 || size < index) larsen_linalg_error ("larsen_linalg_cholesky_insert", "index must be in [0, size].");

	j = index + 1;

	matrix_add_rowcols (size, size, r, 1, 1);

	n = (int) size;
	ldr = n + 1;
	w = (double *) malloc (ldr * sizeof (double));
	dchinx_ (&n, *r, &ldr, &j, u, w, &info);
	free (w);

	return info;
}

void
larsen_linalg_cholesky_delete (const size_t size, double **r, const int index)
{
	int		n;
	int		ldr;
	int		j;
	double	*w;

	if (index < 0 || size <= index) larsen_linalg_error ("larsen_linalg_cholesky_delete", "index must be in [0, size).");

	j = index + 1;

	w = (double *) malloc (size * sizeof (double));

	n = (int) size;
	ldr = (int) size;
	dchdex_ (&n, *r, &ldr, &j, w);
	free (w);

	matrix_remove_rowcols (size, size, r, 1, 1);

	return;
}
