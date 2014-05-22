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
	size_t	p = l->lreg->p;
	double	*z = (double *) malloc (n * sizeof (double));
	double	*y = (double *) malloc (p * sizeof (double));
	/* y(j) = ya(j) for j in A, else = 0 for j not in A */
	for (j = 0; j < p; j++) y[j] = 0.;
	for (j = 0; j < l->sizeA; j++) y[l->A[j]] = ya[j];
	/* z = X * y = X(:, A) * z(A) + X(:, Ac) * 0 */
	dgemv_ ("N", LINREG_CINTP (n), LINREG_CINTP (p), &alpha, x, LINREG_CINTP (n), y, &ione, &dzero, z, &ione);
	free (y);

	return z;
}

/* y = alpha * X(:,A)' * y */
double *
larsen_xa_transpose_dot_y (larsen *l, const size_t n, const double alpha, const double *x, const double *y)
{
	int		j;
	double	*z = (double *) malloc (l->sizeA * sizeof (double));

	/* Evaluate X(:,A)' * y */
	/* another version with dgemv_
	 *
	 *   double	*yp = (double *) malloc (l->p * sizeof (double));
	 *   dgemv_ ("T", CINTP (l->n), CINTP (l->p), &alpha, l->x, CINTP (l->n), z, &ione, &dzero, yp, &ione);
	 *   for (j = 0; j < l->sizeA; j++) y[j] = yp[l->A[j]];
	 *   free (yp);
	 */
	/* The following is faster when l->sizeA is not huge */
	for (j = 0; j < l->sizeA; j++) {
		const double	*xaj = x + LINREG_INDEX_OF_MATRIX (0, l->A[j], n);
		/* z[j] = alpha * X(:,A[j])' * y */
		z[j] = alpha * ddot_ (LINREG_CINTP (n), xaj, &ione, y, &ione);
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
			dcopy_ (LINREG_CINTP (m), *a + LINREG_INDEX_OF_MATRIX (0, j, m), &ione, col, &ione);
			dcopy_ (LINREG_CINTP (m), col, &ione, *a + LINREG_INDEX_OF_MATRIX (0, j, m + dm), &ione);
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
			dcopy_ (LINREG_CINTP (mm), *a + LINREG_INDEX_OF_MATRIX (0, j, m), &ione, col, &ione);
			dcopy_ (LINREG_CINTP (mm), col, &ione, *a + LINREG_INDEX_OF_MATRIX (0, j, m - dm), &ione);
		}
		free (col);
	}

	*a = (double *) realloc (*a, (m - dm) * (n - dn) * sizeof (double));

	return;
}

/* Solve a system of linear equations: L'*L * x = b, where L is a Cholesky factorization. */
int
larsen_linalg_cholesky_svx (const size_t size, double *l, const size_t lda, double *b)
{
	int		info;

	if (!l) linreg_error ("larsen_linalg_cholesky_svx", "first matrix is empty.", __FILE__, __LINE__);
	if (!b) linreg_error ("larsen_linalg_cholesky_svx", "second matrix is empty.", __FILE__, __LINE__);

	dpotrs_ ("U", LINREG_CINTP (size), &ione, l, LINREG_CINTP (lda), b, LINREG_CINTP (size), &info);

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
		linreg_error ("larsen_linalg_cholesky_insert", "index must be in [0, size].", __FILE__, __LINE__);

	j = index + 1;

	matrix_add_rowcols (size, size, r, 1, 1);

	ldr = (int) size + 1;
	w = (double *) malloc (ldr * sizeof (double));
	dchinx_ (LINREG_CINTP (size), *r, &ldr, &j, u, w, &info);
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
		linreg_error ("larsen_linalg_cholesky_delete", "index must be in [0, size).", __FILE__, __LINE__);

	j = index + 1;

	w = (double *) malloc (size * sizeof (double));
	dchdex_ (LINREG_CINTP (size), *r, LINREG_CINTP (size), &j, w);
	free (w);

	matrix_remove_rowcols (size, size, r, 1, 1);

	return;
}

/* Rsolve for economy mode: solve R * x = Q^T * y */
int
larsen_linalg_QR_Rsolve (const size_t size, double *r, double *qty)
{
	int		info;
	if (r == NULL) linreg_error ("c_linalg_QR_Rsolve", "matrix is empty.", __FILE__, __LINE__);
	if (qty == NULL) linreg_error ("c_linalg_QR_Rsolve", "vector is empty.", __FILE__, __LINE__);
	dtrtrs_ ("U", "N", "N", LINREG_CINTP (size), &ione, r, LINREG_CINTP (size), qty, LINREG_CINTP (size), &info);
	return info;
}

/* RTsolve for economy mode: solve R^T * (Q^T * x) = y */
int
larsen_linalg_QR_RTsolve (const size_t size, double *r, double *y)
{
	int		info;
	if (r == NULL) linreg_error ("c_linalg_QR_RTsolve", "matrix is empty.", __FILE__, __LINE__);
	if (y == NULL) linreg_error ("c_linalg_QR_RTsolve", "vector is empty.", __FILE__, __LINE__);
	dtrtrs_ ("U", "T", "N", LINREG_CINTP (size), &ione, r, LINREG_CINTP (size), y, LINREG_CINTP (size), &info);
	return info;
}

/* QR insert for economy mode */
void
larsen_linalg_QR_colinsert (const size_t m, const size_t n, double **q, double **r, const int index, const double *u)
{
	int			j;
	int			k;
	int			ldr;
	int			ldq;
	double		*w;

	if (u == NULL) linreg_error ("larsen_linalg_QR_colinsert", "vector *u is empty.", __FILE__, __LINE__);
	if (index < 0 || n < index) linreg_error ("larsen_linalg_QR_colinsert", "index out of range.", __FILE__, __LINE__);

	k = (m <= n) ? (int) m : (int) n;
	ldq = (int) m;
	ldr = (m <= n + 1) ? (int) m : (int) (n + 1);
	if (m > n) {
		/*-
		 *- 	| Q11 Q12 |    | Q11 Q12 D0 |
		 *- 	| Q21 Q22 | -> | Q21 Q22 D0 |
		 *- 	| Q31 Q32 |    | Q31 Q32 D0 |
		 *-
		 *- 	| R11 R12 |    | R11 R12 D0 |
		 *- 	| D0  R22 | -> | D0  R22 D0 |
		 *- 	               | D0  D0  D0 |
		 */
		matrix_add_rowcols (m, n, q, 0, 1);
		matrix_add_rowcols (n, n, r, 1, 1);
	} else matrix_add_rowcols (m, n, r, 0, 1);

	j = index + 1;	// differ of fortran and C

	w = (double *) malloc (k * sizeof (double));
	dqrinc_ (LINREG_CINTP (m), LINREG_CINTP (n), &k, *q, &ldq, *r, &ldr, &j, u, w);
	free (w);

	return;
}

/* QR delete for economy mode */
void
larsen_linalg_QR_coldelete (const size_t m, const size_t n, double **q, double **r, const int index)
{
	int			j;
	int			k;
	int			ldr;
	int			ldq;
	double		*w;

	if (index < 0 || n <= index) linreg_error ("larsen_linalg_QR_coldelete", "index out of range.", __FILE__, __LINE__);


	k = (m <= n) ? (int) m : (int) n;
	ldq = (int) m;
	ldr = k;

	j = index + 1;	// differ of fortran and C

	w = (double *) malloc ((k - j) * sizeof (double));
	dqrdec_ (LINREG_CINTP (m), LINREG_CINTP (n), &k, *q, &ldq, *r, &ldr, &j, w);
	free (w);

	if (m >= n) {
		/*
		 *- | Q11 Q12 Q13 |    | Q11 Q12 |
		 *- | Q21 Q22 Q23 | -> | Q21 Q22 |
		 *- | Q31 Q32 Q33 |    | Q31 Q32 |
		 *-
		 *- | R11 R12 R13 |    | R11 R12 |
		 *- | D0  R22 R23 | -> | D0  R22 |
		 *- | D0  D0  R33 |
		 */
		matrix_remove_rowcols (m, n, q, 0, 1);
		matrix_remove_rowcols (n, n, r, 1, 1);
	} else matrix_remove_rowcols (m, n, r, 0, 1);

	return;
}
