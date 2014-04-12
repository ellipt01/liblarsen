/*
 * clinalg_cholesky.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <clinalg.h>

/* lapack: cholesky decomposition */
extern void	dpotrf_ (char *uplo, long *n, double *a, long *lda, long *info);
extern void	dpotrs_ (char *uplo, long *n, long *nrhs, double *a, long *lda, double *b, long *ldb, long *info);
/* qrupdate: cholinsert/delete */
extern void	dchinx_ (int *n, double *R, int *ldr, int *j, double *u, double *w, int *info);
extern void	dchdex_ (int *n, double *R, int *ldr, int *j, double *w);

void
clinalg_error (const char * function_name, const char *error_msg)
{
	fprintf (stderr, "ERROR: %s: %s\n", function_name, error_msg);
	exit (1);
}

int
clinalg_cholesky_decomp (size_t size, double *a, size_t lda)
{
	long	info;
	char	uplo;
	long	n;
	long	_lda;

	if (!a) clinalg_error ("clinalg_lapack_dpotrf", "matrix is empty.");

	uplo = 'U';
	n = (long) size;
	_lda = (long) lda;
	dpotrf_ (&uplo, &n, a, &_lda, &info);
	return info;
}

int
clinalg_cholesky_svx (size_t size, double *a, size_t lda, double *b)
{
	long		info;
	char		uplo;
	long		n;
	long		nrhs;
	long		_lda;

	if (!a) clinalg_error ("clinalg_lapack_dpotrs", "first matrix is empty.");
	if (!b) clinalg_error ("clinalg_lapack_dpotrs", "second matrix is empty.");

	uplo = 'U';
	n = (long) size;
	nrhs = 1;
	_lda = (long) lda;
	dpotrs_ (&uplo, &n, &nrhs, a, &_lda, b, &n, &info);

	return info;
}

static void
matrix_add_row_col (size_t m, size_t n, double **a)
{
	int			j;
	double		*col;
	if (a == NULL) {
		*a = (double *) malloc (sizeof (double));
		return;
	}

	*a = (double *) realloc (*a, (m + 1) * (n + 1) * sizeof (double));

	col = (double *) malloc (m * sizeof (double));
	for (j = n; 0 < j; j--) {
		cblas_dcopy (m, *a + index_of_matrix (0, j, m), 1, col, 1);
		cblas_dcopy (m, col, 1, *a + index_of_matrix (0, j, m + 1), 1);
	}
	free (col);

	return;
}

int
clinalg_cholesky_insert (size_t size, double **r, const int index, double *u)
{
	int			info;
	int			n;
	int			ldr;
	int			j;
	double		*w;

	if (index < 0 || size < index) clinalg_error ("clinalg_cholesky_insert", "index must be in [0, size].");

	j = index + 1;

	n = size;
	matrix_add_row_col (n, n, r);

	ldr = n + 1;
	w = (double *) malloc (ldr * sizeof (double));
	dchinx_ (&n, *r, &ldr, &j, u, w, &info);
	free (w);

	return info;
}

static void
matrix_remove_row_col (size_t m, size_t n, double **a)
{
	int			j;
	double		*col;
	if (m <= 1 || n <= 1) {
		if (*a) free (*a);
		return;
	}
	col = (double *) malloc ((m - 1) * sizeof (double));
	for (j = 1; j < n - 1; j++) {
		cblas_dcopy (m - 1, *a + index_of_matrix (0, j, m), 1, col, 1);
		cblas_dcopy (m - 1, col, 1, *a + index_of_matrix (0, j, m - 1), 1);
	}
	free (col);

	*a = (double *) realloc (*a, (m - 1) * (n - 1) * sizeof (double));

	return;
}

void
clinalg_cholesky_delete (size_t size, double **r, const int index)
{
	int		n;
	int		ldr;
	int		j;
	double	*w;

	if (index < 0 || size <= index) clinalg_error ("clinalg_cholesky_delete", "index must be in [0, size).");

	j = index + 1;

	n = size;
	ldr = size;
	w = (double *) malloc (ldr * sizeof (double));

	dchdex_ (&n, *r, &ldr, &j, w);
	free (w);

	matrix_remove_row_col (n, n, r);

	return;
}
