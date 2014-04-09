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

int		c_linalg_lapack_dpotrf (char uplo, c_matrix *a);
int		c_linalg_lapack_dpotrs (char uplo, c_matrix *a, c_matrix *b);

int		c_linalg_cholesky_decomp (c_matrix *a);
int		c_linalg_cholesky_svx (c_matrix *a, c_vector *b);

int		c_linalg_cholesky_insert (c_matrix *r, const int index, const c_vector *u);
void	c_linalg_cholesky_delete (c_matrix *r, const int index);

#ifdef __cplusplus
}
#endif

#endif /* C_LINALG_CHOLESKY_H_ */
