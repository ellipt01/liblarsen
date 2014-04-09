/*
 * c_matrix.h
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#ifndef C_MATRIX_H_
#define C_MATRIX_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define DBL_EPSILON		2.2204460492503131e-16

#define GET_INDEX_OF_VECTOR(v, i) (i * v->stride)
#define GET_INDEX_OF_MATRIX(a, i, j) (i + j * a->lda)

typedef struct s_c_matrix	c_matrix;

struct s_c_matrix {
	size_t		size1;
	size_t		size2;
	size_t		lda;
	size_t		tsize;
	bool		owner;
	double		*data;
};

typedef struct s_c_vector	c_vector;

struct s_c_vector {
	size_t		size;
	size_t		stride;
	size_t		tsize;
	bool		owner;
	double		*data;
};

typedef struct s_c_vector_int	c_vector_int;

struct s_c_vector_int {
	size_t		size;
	size_t		stride;
	size_t		tsize;
	bool		owner;
	int			*data;
};

c_matrix		*c_matrix_alloc (const size_t size1, const size_t size2);
c_vector		*c_vector_alloc (const size_t size);
c_vector_int	*c_vector_int_alloc (const size_t size);

c_matrix		*c_matrix_view_array (const size_t size1, const size_t size2, const size_t lda, double *data);
c_vector		*c_vector_view_array (const size_t size, const size_t stride, double *data);

bool			c_matrix_is_empty (const c_matrix *a);
bool			c_vector_is_empty (const c_vector *a);
bool			c_vector_int_is_empty (const c_vector_int *v);

bool			c_matrix_is_square (const c_matrix *a);

void			c_matrix_free (c_matrix *a);
void			c_vector_free (c_vector *v);
void			c_vector_int_free (c_vector_int *v);

void			c_matrix_set (c_matrix *a, const int i, const int j, double val);
double			c_matrix_get (const c_matrix *a, const int i, const int j);
void			c_vector_set (c_vector *v, const int i, double val);
double			c_vector_get (const c_vector *v, const int i);
void			c_vector_int_set (c_vector_int *v, const int i, int val);
int				c_vector_int_get (const c_vector_int *v, const int i);

void			c_matrix_get_col (c_vector *v, const c_matrix *a, const size_t index);
void			c_matrix_set_col (c_matrix *a, const size_t index, const c_vector *v);

void			c_matrix_add_col (c_matrix *a);
void			c_matrix_add_row (c_matrix *a);
void			c_matrix_remove_col (c_matrix *a);
void			c_matrix_remove_row (c_matrix *a);

void			c_matrix_memcpy (c_matrix *dest, const c_matrix *src);
void			c_vector_memcpy (c_vector *dest, const c_vector *src);

void			c_matrix_set_zero (c_matrix *a);
void			c_vector_set_zero (c_vector *v);

void			c_matrix_fprintf (FILE *stream, const c_matrix *a, const char *format);
void			c_vector_fprintf (FILE *stream, const c_vector *v, const char *format);

/* c_matrixops.c */
void			c_vector_add_constant (c_vector *v, const double x);
double			c_vector_mean (const c_vector *v);
void			c_vector_sub (c_vector *y, const c_vector *x);
double			c_vector_asum (const c_vector *v);
int				c_vector_amax (const c_vector *v);
void			c_vector_scale (c_vector *v, const double alpha);
double			c_vector_nrm (const c_vector *v);
void			c_vector_axpy (double alpha, const c_vector *x, c_vector *y);

double			c_vector_dot_vector (const c_vector *v1, const c_vector *v2);
c_vector		*c_matrix_dot_vector (double alpha, const c_matrix *a, const c_vector *v, double beta);
c_vector		*c_matrix_transpose_dot_vector (double alpha, const c_matrix *a, const c_vector *x, double beta);

#ifdef __cplusplus
}
#endif

#endif /* C_MATRIX_H_ */
