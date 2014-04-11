/*
 * c_linalg.h
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#ifndef C_LINALG_CHOLESKY_H_
#define C_LINALG_CHOLESKY_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <c_matrix.h>

//int		c_linalg_cholesky_decomp (c_matrix *a);
int		c_linalg_cholesky_decomp (size_t size, double *a, size_t lda);
//int		c_linalg_cholesky_svx (c_matrix *a, c_vector *b);
int		c_linalg_cholesky_svx (size_t size, double *a, size_t lda, double *b);


//int		c_linalg_cholesky_insert (c_matrix *r, const int index, double *u);
int		c_linalg_cholesky_insert (size_t n, double **r, const int index, double *u);
//void	c_linalg_cholesky_delete (c_matrix *r, const int index);
void	c_linalg_cholesky_delete (size_t size, double **r, const int index);

#ifdef __cplusplus
}
#endif

#endif /* C_LINALG_CHOLESKY_H_ */
