/*
 * linalg.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <stdlib.h>
#include <larsen.h>

#include "larsen_private.h"

/* y(A(i)) = alpha * xa(i) + y(A(i)) */
void
larsen_axapy (larsen *l, double alpha, double *xa, double *y)
{
	int		i;
	for (i = 0; i < l->sizeA; i++) {
		int		j = l->A[i];
		y[j] += alpha * xa[i];
	}
	return;
}

/* z = alpha * X(:,A) * ya */
double *
larsen_xa_dot_ya (larsen *l, const size_t n, double alpha, const double *x, const double *ya)
{
	int		j;
	size_t	p = l->lsys->p;
	double	*z = (double *) malloc (n * sizeof (double));
	double	*y = (double *) malloc (p * sizeof (double));
	/* y(j) = ya(j) for j in A, else = 0 for j not in A */
	for (j = 0; j < p; j++) y[j] = 0.;
	for (j = 0; j < l->sizeA; j++) y[l->A[j]] = ya[j];
	/* z = X * y = X(:, A) * z(A) + X(:, Ac) * 0 */
	dgemv_ ("N", LINSYS_CINTP (n), LINSYS_CINTP (p), &alpha, x, LINSYS_CINTP (n), y, &ione, &dzero, z, &ione);
	free (y);

	return z;
}

/* y = alpha * X(:,A)' * y */
double *
larsen_xa_transpose_dot_y (larsen *l, const size_t n, const double alpha, const double *x, const double *y)
{
	int		j;
	double	*z = (double *) malloc (l->sizeA * sizeof (double));

	/* eval X(:,A)^T * y */
	/* another version with dgemv_
	 *
	 *   double	*yp = (double *) malloc (l->p * sizeof (double));
	 *   dgemv_ ("T", CINTP (l->n), CINTP (l->p), &alpha, l->x, CINTP (l->n), z, &ione, &dzero, yp, &ione);
	 *   for (j = 0; j < l->sizeA; j++) y[j] = yp[l->A[j]];
	 *   free (yp);
	 *
	 */
	/* The following is faster when l->sizeA is not huge */
	for (j = 0; j < l->sizeA; j++) {
		const double	*xaj = x + LINSYS_INDEX_OF_MATRIX (0, l->A[j], n);
		/* z[j] = alpha * X(:,A[j])' * y */
		z[j] = alpha * ddot_ (LINSYS_CINTP (n), xaj, &ione, y, &ione);
	}
	return z;
}

/* Reallocate *a and add new rows / columns to the end of a matrix.
 * New rows / columns are not initialized.
 * m	:	number of rows of matrix *a
 * n	:	number of columns of matrix *a
 * **a	:	matrix (pointer of m x n array)
 * dm	:	number of rows which will be added (> 0)
 * dn	:	number of columns which will be added (> 0)
 */
static void
matrix_add_rowcols (const size_t m, const size_t n, double **a, const size_t dm, const size_t dn)
{
	int		j;
	if (dm <= 0 && dn <= 0) return;
	if (*a == NULL) {
		*a = (double *) malloc ((m + dm) * (n + dn) * sizeof (double));
		return;
	}

	*a = (double *) realloc (*a, (m + dm) * (n + dn) * sizeof (double));

	if (dm > 0) {
		double	*col = (double *) malloc (m * sizeof (double));
		for (j = n; 0 < j; j--) {
			dcopy_ (LINSYS_CINTP (m), *a + LINSYS_INDEX_OF_MATRIX (0, j, m), &ione, col, &ione);
			dcopy_ (LINSYS_CINTP (m), col, &ione, *a + LINSYS_INDEX_OF_MATRIX (0, j, m + dm), &ione);
		}
		free (col);
	}
	return;
}

/* Remove rows / columns from the end of a matrix and reallocate *a
 * m	:	number of rows of matrix *a
 * n	:	number of columns of matrix *a
 * **a	:	matrix (pointer of m x n array)
 * dm	:	number of rows which will be removed (0 < dm <= m)
 * dn	:	number of columns which will removed (0 < dn <= n)
 */
static void
matrix_remove_rowcols (const size_t m, const size_t n, double **a, const size_t dm, const size_t dn)
{
	int		j;
	if (dm <= 0 && dn <= 0) return;
	if (m - dm <= 0 || n - dn <= 0) {
		if (*a) free (*a);
		*a = NULL;
		return;
	}
	if (dm > 0) {
		double	*col = (double *) malloc ((m - dm) * sizeof (double));
		for (j = 1; j < n - dn; j++) {
			int		mm = (int) (m - dm);
			dcopy_ (LINSYS_CINTP (mm), *a + LINSYS_INDEX_OF_MATRIX (0, j, m), &ione, col, &ione);
			dcopy_ (LINSYS_CINTP (mm), col, &ione, *a + LINSYS_INDEX_OF_MATRIX (0, j, m - dm), &ione);
		}
		free (col);
	}

	*a = (double *) realloc (*a, (m - dm) * (n - dn) * sizeof (double));

	return;
}

/* Solve a system of linear equations: L'*L * x = b, where L is a Cholesky factrization. */
int
larsen_linalg_cholesky_svx (const size_t size, double *l, const size_t lda, double *b)
{
	int		info;

	if (!l) linsys_error ("larsen_linalg_cholesky_svx", "first matrix is empty.", __FILE__, __LINE__);
	if (!b) linsys_error ("larsen_linalg_cholesky_svx", "second matrix is empty.", __FILE__, __LINE__);

	dpotrs_ ("U", LINSYS_CINTP (size), &ione, l, LINSYS_CINTP (lda), b, LINSYS_CINTP (size), &info);

	return info;
}

/* cholinsert: update L -> L1, where L' * L = A and
 * L1' * L1 = [A, u(1:n-1); u(1:n-1)', u(n)] */
int
larsen_linalg_cholesky_insert (const size_t size, double **r, const int index, double *u)
{
	int			info;
	int			ldr;
	int			j;
	double		*w;

	if (index < 0 || size < index)
		linsys_error ("larsen_linalg_cholesky_insert", "index must be in [0, size].", __FILE__, __LINE__);

	j = index + 1;

	matrix_add_rowcols (size, size, r, 1, 1);

	ldr = (int) size + 1;
	w = (double *) malloc (ldr * sizeof (double));
	dchinx_ (LINSYS_CINTP (size), *r, &ldr, &j, u, w, &info);
	free (w);

	return info;
}

/* choldelete: update L -> L1, where
 * L1'*L1 = [
 * 		A(1:index-1,1:index-1),   A(1:index-1,index+1:end);
 * 		A(index+1:end,1:index-1), A(index+1:end, index+1:end)
 * ] */
void
larsen_linalg_cholesky_delete (const size_t size, double **r, const int index)
{
	int		j;
	double	*w;

	if (index < 0 || size <= index)
		linsys_error ("larsen_linalg_cholesky_delete", "index must be in [0, size).", __FILE__, __LINE__);

	j = index + 1;

	w = (double *) malloc (size * sizeof (double));
	dchdex_ (LINSYS_CINTP (size), *r, LINSYS_CINTP (size), &j, w);
	free (w);

	matrix_remove_rowcols (size, size, r, 1, 1);

	return;
}
