/*
 * c_matrix_linlag.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <cblas.h>
#include "c_matrix.h"

extern void	cl_error (const char * function_name, const char *error_msg);
extern int		get_index_of_vector (const c_vector *v, int i);

void
c_vector_add_constant (c_vector *v, const double x)
{
	int		i;
	if (c_vector_is_empty (v)) cl_error ("c_vector_add_constant", "vector *v is empty");
	for (i = 0; i < v->size; i++) v->data[get_index_of_vector (v, i)] += x;
	return;
}

double
c_vector_mean (const c_vector *v)
{
	int	i;
	double	sum;
	if (c_vector_is_empty (v)) cl_error ("c_vector_mean", "vector *v is empty");
	/* x = sum x / N */
	for (i = 0, sum = 0.0; i < v->size; i++) sum += v->data[get_index_of_vector (v, i)];
	return sum / (double) v->size;
}

void
c_vector_sub (c_vector *y, const c_vector *x)
{
	if (c_vector_is_empty (y)) cl_error ("c_vector_sub", "first vector is empty");
	if (c_vector_is_empty (x)) cl_error ("c_vector_sub", "second vector is empty");
	if (x->size != y->size) cl_error ("c_vector_sub", "vector size are not match");
	/* y = y - x */
	cblas_daxpy (y->size, -1.0, x->data, x->stride, y->data, y->stride);
	return;
}

double
c_vector_asum (const c_vector *v)
{
	double asum;
	if (c_vector_is_empty (v)) cl_error ("c_vector_asum", "vector is empty");
	/* x = sum |x| */
	asum = cblas_dasum (v->size, v->data, v->stride);
	return asum;
}

int
c_vector_amax (const c_vector *v)
{
	if (c_vector_is_empty (v)) cl_error ("c_vector_amax", "vector is empty");
	return cblas_idamax (v->size, v->data, v->stride);
}

void
c_vector_scale (c_vector *v, const double alpha)
{
	if (c_vector_is_empty (v)) cl_error ("c_vector_scale", "vector is empty");
	/* x = alpha * x */
	cblas_dscal (v->size, alpha, v->data, v->stride);
	return;
}

double
c_vector_nrm (const c_vector *v)
{
	if (c_vector_is_empty (v)) cl_error ("c_vector_nrm", "vector is empty");
	return cblas_dnrm2 (v->size, v->data, v->stride);
}

void
c_vector_axpy (const double alpha, const c_vector *x, c_vector *y)
{
	if (c_vector_is_empty (x)) cl_error ("c_vector_axpy", "first vector is empty");
	if (c_vector_is_empty (y)) cl_error ("c_vector_axpy", "second vector is empty");
	if (x->size != y->size) cl_error ("c_vector_axpy", "vector size are not match");
	/* y = y + alpha * x */
	cblas_daxpy (y->size, alpha, x->data, x->stride, y->data, y->stride);
	return;
}

double
c_vector_dot_vector (const c_vector *v1, const c_vector *v2)
{
	if (c_vector_is_empty (v1)) cl_error ("c_vector_dot_vector", "first vector is empty");
	if (c_vector_is_empty (v2)) cl_error ("c_vector_dot_vector", "second vector is empty");
	return cblas_ddot (v1->size, v1->data, v1->stride, v2->data, v2->stride);
}

c_vector *
c_matrix_dot_vector (const c_matrix *a, const c_vector *v)
{
	c_vector *r;
	if (c_matrix_is_empty (a)) cl_error ("c_matrix_dot_vector", "matrix is empty");
	if (c_vector_is_empty (v)) cl_error ("c_matrix_dot_vector", "vector is empty");
	if (a->size2 != v->size) cl_error ("c_matrix_dot_vector", "vector and matrix size are not match");
	r = c_vector_alloc (a->size1);
	cblas_dgemv (CblasColMajor, CblasNoTrans, a->size1, a->size2, 1.0, a->data, a->lda, v->data, v->stride, 0.0, r->data, r->stride);
	return r;
}

c_vector *
c_matrix_transpose_dot_vector (const c_matrix *a, const c_vector *v)
{
	c_vector *r;
	if (c_matrix_is_empty (a)) cl_error ("c_matrix_transpose_dot_vector", "matrix is empty");
	if (c_vector_is_empty (v)) cl_error ("c_matrix_transpose_dot_vector", "vector is empty");
	if (a->size1 != v->size) cl_error ("c_matrix_transpose_dot_vector", "vector and matrix size are not match");
	r = c_vector_alloc (a->size2);
	cblas_dgemv (CblasColMajor, CblasTrans, a->size1, a->size2, 1.0, a->data, a->lda, v->data, v->stride, 0.0, r->data, r->stride);
	return r;
}

