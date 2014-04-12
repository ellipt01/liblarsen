/*
 * c_linalg.h
 *
 *  Wrapper of lapack, blas and qrupdate with a GSL like interface
 *
 *  Created on: 2014/04/08
 *      Author: utsugi
 */

#ifndef C_LINALG_H_
#define C_LINALG_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <cblas.h>

#define index_of_matrix(i, j, lda) ((i) + (j) * (lda))

int		c_linalg_cholesky_decomp (size_t size, double *a, size_t lda);
int		c_linalg_cholesky_svx (size_t size, double *a, size_t lda, double *b);

int		c_linalg_cholesky_insert (size_t n, double **r, const int index, double *u);
void	c_linalg_cholesky_delete (size_t size, double **r, const int index);

#ifdef __cplusplus
}
#endif

#endif /* C_LINALG_H_ */
