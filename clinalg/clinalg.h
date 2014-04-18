/*
 * clinalg.h
 *
 *  Wrapper of lapack and qrupdate
 *
 *  Created on: 2014/04/08
 *      Author: utsugi
 */

#ifndef CLINALG_H_
#define CLINALG_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <cblas.h>

#ifndef index_of_matrix
#define index_of_matrix(i, j, lda) ((i) + (j) * (lda))
#endif

/* lapack: cholesky decomposition */
#ifndef HAVE_LAPACK_H
extern void	dpotrf_ (char *uplo, int *n, double *a, int *lda, int *info);
extern void	dpotrs_ (char *uplo, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *info);
#endif

/* qrupdate: cholinsert/delete */
#ifndef HAVE_QRUPDATE_H
extern void	dchinx_ (int *n, double *R, int *ldr, int *j, double *u, double *w, int *info);
extern void	dchdex_ (int *n, double *R, int *ldr, int *j, double *w);
#endif

int		clinalg_cholesky_decomp (const size_t size, double *a, const size_t lda);
int		clinalg_cholesky_svx (const size_t size, double *a, const size_t lda, double *b);

int		clinalg_cholesky_insert (const size_t n, double **r, const int index, double *u);
void	clinalg_cholesky_delete (const size_t size, double **r, const int index);

#ifdef __cplusplus
}
#endif

#endif /* CLINALG_H_ */
