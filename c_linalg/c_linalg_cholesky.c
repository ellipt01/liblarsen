/*
 * cl_linalg.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include "c_matrix.h"

extern void	cl_error (const char * function_name, const char *error_msg);

void	dpotrf_ (char *uplo, long *n, double *a, long *lda, long *info);
void	dpotrs_ (char *uplo, long *n, long *nrhs, double *a, long *lda, double *b, long *ldb, long *info);
void	dchinx_ (int *n, double *R, int *ldr, int *j, double *u, double *w, int *info);
void	dchdex_ (int *n, double *R, int *ldr, int *j, double *w);

int
cl_linalg_cholesky_decomp (c_matrix *a)
{
	long		info;
	char		uplo;
	long		n;
	long		lda;

	if (c_matrix_is_empty (a)) cl_error ("cl_linalg_cholesky_decomp", "matrix *a is empty");
	if (!c_matrix_is_square (a)) cl_error ("cl_linalg_cholesky_decomp", "matrix *a must be square");

	uplo = 'U';
	n = (long) a->size2;
	lda = (long) a->lda;
	dpotrf_ (&uplo, &n, a->data, &lda, &info);
	return (int) info;
}

int
cl_linalg_cholesky_svx (c_matrix *a, c_vector *b)
{
	long		info;
	char		uplo;
	long		n;
	long		nrhs;
	long		lda;
	long		ldb;

	if (c_matrix_is_empty (a)) cl_error ("cl_linalg_cholesky_svx", "matrix *a is empty");
	if (c_vector_is_empty (b)) cl_error ("cl_linalg_cholesky_svx", "vector *b is empty");
	if (!c_matrix_is_square (a)) cl_error ("cl_linalg_cholesky_svx", "matrix *a must be square");

	uplo = 'U';
	n = (long) a->size1;
	nrhs = 1;
	lda = (long) a->lda;
	ldb = (long) b->size;
	dpotrs_ (&uplo, &n, &nrhs, a->data, &lda, b->data, &ldb, &info);
	return (int) info;
}

int
cl_linalg_cholesky_insert (c_matrix *r, const int index, const c_vector *u)
{
	int		n;
	int		ldr;
	int		j;
	double	*w;
	int		info;

	if (c_matrix_is_empty (r)) cl_error ("cl_linalg_cholesky_insert", "matrix *r is empty");
	if (c_vector_is_empty (u)) cl_error ("cl_linalg_cholesky_insert", "vector *u is empty");
	if (r->size1 != r->size2) cl_error ("cl_linalg_cholesky_insert", "matrix *r must be square");
	if (u->size != r->size1 + 1) cl_error ("cl_linalg_cholesky_insert", "size of matrix *r and vector *u are not match");
	if (index < 0 || r->size1 < index) cl_error ("cl_linalg_cholesky_insert", "index must be in [0, r->size1]");

	j = index + 1;

	n = r->size1;

	c_matrix_add_col (r);
	c_matrix_add_row (r);

	ldr = r->lda;
	w = (double *) malloc (ldr * sizeof (double));
	dchinx_ (&n, r->data, &ldr, &j, u->data, w, &info);
	free (w);

	return info;
}

void
cl_linalg_cholesky_delete (c_matrix *r, const int index)
{
	int		n;
	int		ldr;
	int		j;
	double	*w;

	if (c_matrix_is_empty (r)) cl_error ("cl_linalg_cholesky_delete", "matrix *r is empty");
	if (r->size1 != r->size2) cl_error ("cl_linalg_cholesky_delete", "matrix *r must be square");
	if (index < 0 || r->size1 <= index) cl_error ("cl_linalg_cholesky_delete", "index must be in [0, r->size1)");

	j = index + 1;

	n = r->size1;
	ldr = r->lda;
	w = (double *) malloc (ldr * sizeof (double));

	dchdex_ (&n, r->data, &ldr, &j, w);
	free (w);
	c_matrix_remove_col (r);
	if (!c_matrix_is_empty (r)) c_matrix_remove_row (r);

	return;
}
