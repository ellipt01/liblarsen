/*
 * cl_matrix.h
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#ifndef CL_MATRIX_H_
#define CL_MATRIX_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define DBL_EPSILON		2.2204460492503131e-16

#define CL_MIN(a, b)		((a > b) ? b : a)

typedef struct s_cl_matrix	cl_matrix;

struct s_cl_matrix {
	size_t		size1;
	size_t		size2;
	size_t		lda;
	size_t		tsize;
	bool		owner;
	bool		allocated;
	double		*data;
};

typedef struct s_cl_vector	cl_vector;

struct s_cl_vector {
	size_t		size;
	size_t		stride;
	size_t		tsize;
	bool		owner;
	bool		allocated;
	double		*data;
};

typedef struct s_cl_vector_int	cl_vector_int;

struct s_cl_vector_int {
	size_t		size;
	size_t		stride;
	size_t		tsize;
	bool		owner;
	bool		allocated;
	int			*data;
};

cl_matrix		*cl_matrix_alloc (const size_t size1, const size_t size2);
cl_vector		*cl_vector_alloc (const size_t size);
cl_vector_int	*cl_vector_int_alloc (const size_t size);

bool			cl_matrix_is_empty (const cl_matrix *a);
bool			cl_vector_is_empty (const cl_vector *a);
bool			cl_vector_int_is_empty (const cl_vector_int *v);

bool			cl_matrix_is_square (const cl_matrix *a);

void			cl_matrix_free (cl_matrix *a);
void			cl_vector_free (cl_vector *v);
void			cl_vector_int_free (cl_vector_int *v);

void			cl_matrix_set (cl_matrix *a, const int i, const int j, double val);
double			cl_matrix_get (const cl_matrix *a, const int i, const int j);
void			cl_vector_set (cl_vector *v, const int i, double val);
double			cl_vector_get (const cl_vector *v, const int i);
void			cl_vector_int_set (cl_vector_int *v, const int i, int val);
int				cl_vector_int_get (const cl_vector_int *v, const int i);

void			cl_matrix_get_col (cl_vector *v, const cl_matrix *a, const size_t index);
void			cl_matrix_set_col (cl_matrix *a, const size_t index, const cl_vector *v);

void			cl_matrix_add_col (cl_matrix *a);
void			cl_matrix_add_row (cl_matrix *a);
void			cl_matrix_remove_col (cl_matrix *a);
void			cl_matrix_remove_row (cl_matrix *a);

void			cl_matrix_memcpy (cl_matrix *dest, const cl_matrix *src);
void			cl_vector_memcpy (cl_vector *dest, const cl_vector *src);

void			cl_matrix_set_zero (cl_matrix *a);
void			cl_vector_set_zero (cl_vector *v);

void			cl_matrix_fprintf (FILE *stream, const cl_matrix *a, const char *format);
void			cl_vector_fprintf (FILE *stream, const cl_vector *v, const char *format);

/* cl_matrix_ops.c */
void			cl_vector_add_constant (cl_vector *v, const double x);
double			cl_vector_mean (const cl_vector *v);
void			cl_vector_sub (cl_vector *y, const cl_vector *x);
double			cl_vector_asum (const cl_vector *v);
int				cl_vector_amax (const cl_vector *v);
void			cl_vector_scale (cl_vector *v, const double alpha);
double			cl_vector_nrm (const cl_vector *v);
void			cl_vector_axpy (const double alpha, const cl_vector *x, cl_vector *y);

double			cl_vector_dot_vector (const cl_vector *v1, const cl_vector *v2);
cl_vector		*cl_matrix_dot_vector (const cl_matrix *a, const cl_vector *v);
cl_vector		*cl_matrix_transpose_dot_vector (const cl_matrix *a, const cl_vector *v);

#ifdef __cplusplus
}
#endif

#endif /* CL_MATRIX_H_ */
