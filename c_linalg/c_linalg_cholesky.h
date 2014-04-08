/*
 * cl_linalg.h
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#ifndef CL_LINALG_H_
#define CL_LINALG_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <c_matrix.h>

int		cl_linalg_cholesky_decomp (c_matrix *a);
int		cl_linalg_cholesky_svx (c_matrix *a, c_vector *b);

int		cl_linalg_cholesky_insert (c_matrix *r, const int index, const c_vector *u);
void	cl_linalg_cholesky_delete (c_matrix *r, const int index);

#ifdef __cplusplus
}
#endif

#endif /* CL_LINALG_H_ */
