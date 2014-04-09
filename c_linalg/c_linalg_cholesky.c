/*
 * c_linalg.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include "c_matrix.h"

extern void	cl_error (const char * function_name, const char *error_msg);

/* lapack: cholesky decomposition */
extern void	dpotrf_ (char *uplo, long *n, double *a, long *lda, long *info);
extern void	dpotrs_ (char *uplo, long *n, long *nrhs, double *a, long *lda, double *b, long *ldb, long *info);
/* qrupdate: cholinsert/delete */
extern void	dchinx_ (int *n, double *R, int *ldr, int *j, double *u, double *w, int *info);
extern void	dchdex_ (int *n, double *R, int *ldr, int *j, double *w);

/* interface of lapack dpotrf_ */
int
c_linalg_lapack_dpotrf (char uplo, c_matrix *a)
{
	long	info;
	long	n;
	long	lda;

	if (c_matrix_is_empty (a)) cl_error ("c_linalg_lapack_dpotrf", "matrix is empty");
	if (!c_matrix_is_square (a)) cl_error ("c_linalg_lapack_dpotrf", "matrix must be square");

	n = (long) a->size1;
	lda = (long) a->lda;
	dpotrf_ (&uplo, &n, a->data, &lda, &info);
	return (int) info;
}

/* interface of lapack dpotrs_ */
int
c_linalg_lapack_dpotrs (char uplo, c_matrix *a, c_matrix *b)
{
	long		info;
	long		n;
	long		nrhs;
	long		lda;
	long		ldb;

	if (c_matrix_is_empty (a)) cl_error ("c_linalg_lapack_dpotrs", "first matrix is empty");
	if (c_matrix_is_empty (b)) cl_error ("c_linalg_lapack_dpotrs", "second matrix is empty");
	if (!c_matrix_is_square (a)) cl_error ("c_linalg_lapack_dpotrs", "first matrix must be square");
	if (a->size1 != b->size1) cl_error ("c_linalg_lapack_dpotrs", "matrix size does not match");

	n = (long) a->size1;
	nrhs = (long) b->size2;
	lda = (long) a->lda;
	ldb = (long) b->lda;
	dpotrs_ (&uplo, &n, &nrhs, a->data, &lda, b->data, &ldb, &info);
	return (int) info;
}

int
c_linalg_cholesky_decomp (c_matrix *a)
{
	int			info;
	if (c_matrix_is_empty (a)) cl_error ("c_linalg_cholesky_decomp", "matrix is empty");
	if (!c_matrix_is_square (a)) cl_error ("c_linalg_cholesky_decomp", "matrix must be square");

	info = c_linalg_lapack_dpotrf ('U', a);
	return info;
}

int
c_linalg_cholesky_svx (c_matrix *a, c_vector *b)
{
	int			info;
	c_matrix	*c;
	if (c_matrix_is_empty (a)) cl_error ("c_linalg_cholesky_svx", "first matrix is empty");
	if (c_vector_is_empty (b)) cl_error ("c_linalg_cholesky_svx", "second vector is empty");
	if (!c_matrix_is_square (a)) cl_error ("c_linalg_cholesky_svx", "first matrix must be square");
	if (a->size1 != b->size) cl_error ("c_linalg_cholesky_svx", "matrix size does not match");

	c = c_matrix_view_array (b->size, 1, b->size, b->data);
	info = c_linalg_lapack_dpotrs ('U', a, c);
	c_matrix_free (c);
	return info;
}

int
c_linalg_cholesky_insert (c_matrix *r, const int index, const c_vector *u)
{
	int		n;
	int		ldr;
	int		j;
	double	*w;
	int		info;

	if (c_matrix_is_empty (r)) cl_error ("c_linalg_cholesky_insert", "matrix is empty");
	if (c_vector_is_empty (u)) cl_error ("c_linalg_cholesky_insert", "vector is empty");
	if (r->size1 != r->size2) cl_error ("c_linalg_cholesky_insert", "matrix must be square");
	if (u->size != r->size1 + 1) cl_error ("c_linalg_cholesky_insert", "matrix and vector size does not match");
	if (index < 0 || r->size1 < index) cl_error ("c_linalg_cholesky_insert", "index must be in [0, r->size1]");

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

int
c_linalg_cholesky_delete (c_matrix *r, const int index)
{
	int		info = 0;
	int		n;
	int		ldr;
	int		j;
	double	*w;

	if (c_matrix_is_empty (r)) cl_error ("c_linalg_cholesky_delete", "matrix is empty");
	if (r->size1 != r->size2) cl_error ("c_linalg_cholesky_delete", "matrix must be square");
	if (index < 0 || r->size1 <= index) cl_error ("c_linalg_cholesky_delete", "index must be in [0, r->size1)");

	j = index + 1;

	n = r->size1;
	ldr = r->lda;
	w = (double *) malloc (ldr * sizeof (double));

	dchdex_ (&n, r->data, &ldr, &j, w);
	free (w);
	c_matrix_remove_col (r);
	if (!c_matrix_is_empty (r)) c_matrix_remove_row (r);

	return info;
}
