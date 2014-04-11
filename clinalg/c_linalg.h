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

#include <c_matrix.h>
#include <c_linalg_cholesky.h>

#define index_of_matrix(i, j, lda) ((i) + (j) * (lda))

#ifdef __cplusplus
}
#endif

#endif /* C_LINALG_H_ */
