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
int		c_linalg_cholesky_decomp (size_t size, double *a);
//int		c_linalg_cholesky_svx (c_matrix *a, c_vector *b);
int		c_linalg_cholesky_svx (size_t size, double *a, double *b);

int		c_linalg_cholesky_insert (int n, double *r, const int index, double *u);
void	c_linalg_cholesky_delete (int n, double *r, const int index);

#ifdef __cplusplus
}
#endif

#endif /* C_LINALG_CHOLESKY_H_ */
