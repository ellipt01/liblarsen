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

#define index_of_matrix(i, j, lda) ((i) + (j) * (lda))

int		clinalg_cholesky_decomp (size_t size, double *a, size_t lda);
int		clinalg_cholesky_svx (size_t size, double *a, size_t lda, double *b);

int		clinalg_cholesky_insert (size_t n, double **r, const int index, double *u);
void	clinalg_cholesky_delete (size_t size, double **r, const int index);

#ifdef __cplusplus
}
#endif

#endif /* CLINALG_H_ */
