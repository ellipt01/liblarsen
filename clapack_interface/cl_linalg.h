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

int		cl_linalg_cholesky_decomp (cl_matrix *a);
int		cl_linalg_cholesky_svx (cl_matrix *a, cl_vector *b);

int		cl_linalg_cholesky_insert (cl_matrix *r, const int index, const cl_vector *u);
void	cl_linalg_cholesky_delete (cl_matrix *r, const int index);

#ifdef __cplusplus
}
#endif

#endif /* CL_LINALG_H_ */
