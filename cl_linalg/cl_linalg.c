/*
 * cl_linalg.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include "cl_matrix.h"

extern void	cl_error (const char * function_name, const char *error_msg);

void	dpotrf_ (char *uplo, long *n, double *a, long *lda, long *info);
void	dpotrs_ (char *uplo, long *n, long *nrhs, double *a, long *lda, double *b, long *ldb, long *info);
void	dchinx_ (int *n, double *R, int *ldr, int *j, double *u, double *w, int *info);
void	dchdex_ (int *n, double *R, int *ldr, int *j, double *w);


int
cl_linalg_cholesky_decomp (cl_matrix *a)
{
	long		info;

	if (cl_matrix_is_empty (a)) cl_error ("cl_linalg_cholesky_decomp", "matrix A is empty");
	if (!cl_matrix_is_square (a)) cl_error ("cl_linalg_cholesky_decomp", "matrix A must be square");
	{
		int		i, j;
		char	uplo = 'U';
		long	n = (long) a->size2;
		long	lda = (long) a->lda;
		dpotrf_ (&uplo, &n, a->data, &lda, &info);
		for (i = 0; i < a->size2; i++) {
			for (j = 0; j < i; j++) cl_matrix_set (a, i, j, cl_matrix_get (a, j, i));
		}
	}
	return (int) info;
}

static int
cl_linalg_lapack_dpotrs (char uplo, cl_matrix *a, cl_matrix *b)
{
	long		info;

	long		n;
	long		nrhs;
	long		lda;
	long		ldb;

	if (!cl_matrix_is_square (a)) cl_error ("cl_linalg_lapack_dpotrs", "matrix *a must be square: *a must be cholesky decomposed matrix");
	if (uplo != 'L' && uplo != 'U') cl_error ("cl_linalg_lapack_dpotrs", "uplo must be 'L' or 'U'");
	n = (long) a->size1;
	nrhs = (long) b->size2;
	lda = (long) a->lda;
	ldb = (long) b->lda;
	dpotrs_ (&uplo, &n, &nrhs, a->data, &lda, b->data, &ldb, &info);

	return (int) info;
}

int
cl_linalg_cholesky_svx (cl_matrix *a, cl_vector *b)
{
	int			info;
	cl_matrix	*x;

	if (cl_matrix_is_empty (a)) cl_error ("cl_linalg_cholesky_svx", "matrix *a is empty");
	if (cl_vector_is_empty (b)) cl_error ("cl_linalg_cholesky_svx", "vector *b is empty");
	if (!cl_matrix_is_square (a)) cl_error ("cl_linalg_cholesky_svx", "matrix *a must be square");

	x = cl_matrix_view_array (b->size, 1, b->data);
	info = cl_linalg_lapack_dpotrs ('U', a, x);
	cl_matrix_free (x);
	return info;
}



int
cl_linalg_cholesky_insert (cl_matrix *r, int index, cl_vector *u)
{
	int		n;
	int		ldr;
	int		j;
	double	*w;
	int		info;

	if (cl_matrix_is_empty (r)) cl_error ("cl_linalg_cholesky_insert", "matrix *r is empty");
	if (cl_vector_is_empty (u)) cl_error ("cl_linalg_cholesky_insert", "vector *u is empty");
	if (r->size1 != r->size2) cl_error ("cl_linalg_cholesky_insert", "matrix *r must be square");
	if (u->size != r->size1 + 1) cl_error ("cl_linalg_cholesky_insert", "size of matrix *r and vector *u invalid");
	if (index < 0) cl_error ("cl_linalg_cholesky_insert", "index must be >= 0");
	if (index > r->size1) cl_error ("cl_linalg_cholesky_insert", "index must be <= r->size1");

	j = index + 1;

	n = r->size1;

	cl_matrix_add_col (r);
	cl_matrix_add_row (r);

	ldr = r->lda;
	w = (double *) malloc (ldr * sizeof (double));
	dchinx_ (&n, r->data, &ldr, &j, u->data, w, &info);
	free (w);

	return info;
}

void
cl_linalg_cholesky_delete (cl_matrix *r, int index)
{
	int		n;
	int		ldr;
	int		j;
	double	*w;

	if (cl_matrix_is_empty (r)) cl_error ("cl_linalg_cholesky_delete", "matrix *r is empty");
	if (r->size1 != r->size2) cl_error ("cl_linalg_cholesky_delete", "matrix *r must be square");
	if (index < 0) cl_error ("cl_linalg_cholesky_delete", "index must be >= 0");
	if (index >= r->size1) cl_error ("cl_linalg_cholesky_delete", "index must be < r->size1");

	j = index + 1;

	n = r->size1;
	ldr = r->lda;
	w = (double *) malloc (ldr * sizeof (double));

	dchdex_ (&n, r->data, &ldr, &j, w);
	free (w);
	cl_matrix_remove_col (r);
	if (!cl_matrix_is_empty (r)) cl_matrix_remove_row (r);

	return;
}
