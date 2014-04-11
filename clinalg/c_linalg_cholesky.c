/*
 * c_linalg_cholesky.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <c_linalg.h>

/* c_matrix.c */
extern void	c_error (const char * function_name, const char *error_msg);

/* lapack: cholesky decomposition */
extern void	dpotrf_ (char *uplo, long *n, double *a, long *lda, long *info);
extern void	dpotrs_ (char *uplo, long *n, long *nrhs, double *a, long *lda, double *b, long *ldb, long *info);
/* qrupdate: cholinsert/delete */
extern void	dchinx_ (int *n, double *R, int *ldr, int *j, double *u, double *w, int *info);
extern void	dchdex_ (int *n, double *R, int *ldr, int *j, double *w);

int
c_linalg_cholesky_decomp (size_t size, double *a)
{
	long	info;
	char	uplo;
	long	n;
	long	lda;

	if (!a) c_error ("c_linalg_lapack_dpotrf", "matrix is empty.");

	uplo = 'U';
	n = (long) size;
	lda = (long) size;
	dpotrf_ (&uplo, &n, a, &lda, &info);
	return info;
}

int
c_linalg_cholesky_svx (size_t size, double *a, double *b)
{
	long		info;
	char		uplo;
	long		n;
	long		nrhs;
	long		lda;
	long		ldb;

	if (!a) c_error ("c_linalg_lapack_dpotrs", "first matrix is empty.");
	if (!b) c_error ("c_linalg_lapack_dpotrs", "second matrix is empty.");

	uplo = 'U';
	n = (long) size;
	nrhs = 1;
	lda = (long) size;
	ldb = (long) size;
	dpotrs_ (&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);

	return info;
}

static void
matrix_add_row (size_t size1, size_t size2, double *a)
{
	int			j;
	size_t		m = size1 + 1;
	size_t		n = size2;
	double	*col = (double *) malloc (m * sizeof (double));

	a = (double *) realloc (a, m * n * sizeof (double));

	for (j = n - 1; 0 < j; j--) {
		cblas_dcopy (m, a + j * m, 1, col, 1);
		cblas_dcopy (m, col, 1, a + j * m, 1);
	}
	free (col);
	for (j = 0; j < n; j++) a[index_of_matrix (m - 1, j, m)] = 0;
	return;
}

static void
matrix_add_col (size_t size1, size_t size2, double *a)
{
	size_t		m = size1;
	size_t		n = size2 + 1;
	a = (double *) realloc (a, m * n * sizeof (double));
	return;
}

static void
matrix_remove_row (size_t size1, size_t size2, double *a)
{
	size_t		m = size1 - 1;
	size_t		n = size2;

	if (m > 0) {
		int		j;
		double	*col = (double *) malloc (m * sizeof (double));
		for (j = 1; j < n; j++) {
			cblas_dcopy (n, a + j * m, 1, col, 1);
			cblas_dcopy (n, col, 1, a + j * m, 1);
		}
		free (col);
	}
	if (a) a = (double *) realloc (a, m * n * sizeof (double));

	return;
}

static void
matrix_remove_col (size_t size1, size_t size2, double *a)
{
	size_t		m = size1;
	size_t		n = size2 - 1;

	if (a) a = (double *) realloc (a, m * n * sizeof (double));

	return;
}

int
c_linalg_cholesky_insert (int n, double *r, const int index, double *u)
{
	int		j;
	double	*w;
	int		info;

	if (!r) c_error ("c_linalg_cholesky_insert", "matrix is empty.");
	if (!u) c_error ("c_linalg_cholesky_insert", "vector is empty.");
	if (index < 0 || n < index) c_error ("c_linalg_cholesky_insert", "index must be in [0, n].");

	j = index + 1;

	matrix_add_col (n, n, r);
	matrix_add_row (n + 1, n, r);

	w = (double *) malloc (n * sizeof (double));
	dchinx_ (&n, r, &n, &j, u, w, &info);
	free (w);

	return info;
}

void
c_linalg_cholesky_delete (int n, double *r, const int index)
{
	int		j;
	double	*w;

	if (!r) c_error ("c_linalg_cholesky_delete", "matrix is empty.");
	if (index < 0 || n <= index) c_error ("c_linalg_cholesky_delete", "index must be in [0, n).");

	j = index + 1;

	w = (double *) malloc (n * sizeof (double));

	dchdex_ (&n, r, &n, &j, w);
	free (w);
	matrix_remove_col (n, n, r);
	if (r) matrix_remove_row (n, n - 1, r);

	return;
}
