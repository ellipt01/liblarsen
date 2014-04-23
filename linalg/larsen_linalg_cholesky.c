/*
 * larsen_linalg_cholesky.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <larsen_linalg.h>

void
larsen_linalg_error (const char * function_name, const char *error_msg)
{
	fprintf (stderr, "ERROR: %s: %s\n", function_name, error_msg);
	exit (1);
}

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

static void
matrix_add_row_col (const size_t m, const size_t n, double **a)
{
	int			j;
	double		*col;
	if (m == 0 || n == 0) {
		*a = (double *) malloc (sizeof (double));
		return;
	}

	*a = (double *) realloc (*a, (m + 1) * (n + 1) * sizeof (double));

	col = (double *) malloc (m * sizeof (double));
	for (j = n; 0 < j; j--) {
		cblas_dcopy (m, *a + INDEX_OF_MATRIX (0, j, m), 1, col, 1);
		cblas_dcopy (m, col, 1, *a + INDEX_OF_MATRIX (0, j, m + 1), 1);
	}
	free (col);

	return;
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

	matrix_add_row_col (size, size, r);

	n = (int) size;
	ldr = n + 1;
	w = (double *) malloc (ldr * sizeof (double));
	dchinx_ (&n, *r, &ldr, &j, u, w, &info);
	free (w);

	return info;
}

static void
matrix_remove_row_col (const size_t m, const size_t n, double **a)
{
	int			j;
	double		*col;
	if (m <= 1 || n <= 1) {
		if (*a) free (*a);
		*a = NULL;
		return;
	}
	col = (double *) malloc ((m - 1) * sizeof (double));
	for (j = 1; j < n - 1; j++) {
		cblas_dcopy (m - 1, *a + INDEX_OF_MATRIX (0, j, m), 1, col, 1);
		cblas_dcopy (m - 1, col, 1, *a + INDEX_OF_MATRIX (0, j, m - 1), 1);
	}
	free (col);

	*a = (double *) realloc (*a, (m - 1) * (n - 1) * sizeof (double));

	return;
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

	matrix_remove_row_col (size, size, r);

	return;
}
