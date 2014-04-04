/*
 * cl_matrix_linlag.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <cblas.h>
#include "cl_matrix.h"

extern void	cl_error (const char * function_name, const char *error_msg);
extern int		get_index_of_vector (const cl_vector *v, int i);

void
cl_vector_add_constant (cl_vector *v, const double x)
{
	int		i;
	if (cl_vector_is_empty (v)) cl_error ("cl_vector_add_constant", "vector v is empty");
	for (i = 0; i < v->size; i++) v->data[get_index_of_vector (v, i)] += x;
	return;
}

double
cl_vector_mean (const cl_vector *v)
{
	int	i;
	double	sum;
	if (cl_vector_is_empty (v)) cl_error ("cl_vector_mean", "vector is empty");
	/* x = sum x / N */
	for (i = 0, sum = 0.0; i < v->size; i++) sum += v->data[get_index_of_vector (v, i)];
	return sum / (double) v->size;
}

void
cl_vector_sub (cl_vector *y, const cl_vector *x)
{
	if (cl_vector_is_empty (y)) cl_error ("cl_vector_sub", "first vector is empty");
	if (cl_vector_is_empty (x)) cl_error ("cl_vector_sub", "second vector is empty");
	if (x->size != y->size) cl_error ("cl_vector_sub", "vector size not match");
	/* y = y - x */
	cblas_daxpy (y->size, -1.0, x->data, x->stride, y->data, y->stride);
	return;
}

double
cl_vector_asum (const cl_vector *v)
{
	double asum;
	if (cl_vector_is_empty (v)) cl_error ("cl_vector_asum", "vector is empty");
	/* x = sum |x| */
	asum = cblas_dasum (v->size, v->data, v->stride);
	return asum;
}

int
cl_vector_amax (const cl_vector *v)
{
	if (cl_vector_is_empty (v)) cl_error ("cl_vector_amax", "vector is empty");
	return cblas_idamax (v->size, v->data, v->stride);
}

void
cl_vector_scale (cl_vector *v, const double alpha)
{
	if (cl_vector_is_empty (v)) cl_error ("cl_vector_scale", "vector is empty");
	/* x = alpha * x */
	cblas_dscal (v->size, alpha, v->data, v->stride);
	return;
}

double
cl_vector_nrm (const cl_vector *v)
{
	if (cl_vector_is_empty (v)) cl_error ("cl_vector_nrm", "vector is empty");
	return cblas_dnrm2 (v->size, v->data, v->stride);
}

void
cl_vector_axpy (const double alpha, const cl_vector *x, cl_vector *y)
{
	if (cl_vector_is_empty (x)) cl_error ("cl_vector_axpy", "first vector is empty");
	if (cl_vector_is_empty (y)) cl_error ("cl_vector_axpy", "second vector is empty");
	if (x->size != y->size) cl_error ("cl_vector_axpy", "vector size not match");
	/* y = y + alpha * x */
	cblas_daxpy (y->size, alpha, x->data, x->stride, y->data, y->stride);
	return;
}

double
cl_vector_dot_vector (const cl_vector *v1, const cl_vector *v2)
{
	if (cl_vector_is_empty (v1)) cl_error ("cl_vector_dot_vector", "vector v1 is empty");
	if (cl_vector_is_empty (v2)) cl_error ("cl_vector_dot_vector", "vector v2 is empty");
	return cblas_ddot (v1->size, v1->data, v1->stride, v2->data, v2->stride);
}

cl_vector *
cl_matrix_dot_vector (const cl_matrix *a, const cl_vector *v)
{
	cl_vector *r;
	if (cl_matrix_is_empty (a)) cl_error ("cl_matrix_dot_vector", "matrix is empty");
	if (cl_vector_is_empty (v)) cl_error ("cl_matrix_dot_vector", "vector is empty");
	if (a->size2 != v->size) cl_error ("cl_matrix_dot_vector", "vector and matrix size not match");
	r = cl_vector_alloc (a->size1);
	if (cl_vector_is_empty (r)) cl_error ("cl_matrix_dot_vector", "vector not allocated");
	cblas_dgemv (CblasColMajor, CblasNoTrans, a->size1, a->size2, 1.0, a->data, a->lda, v->data, v->stride, 0.0, r->data, r->stride);
	return r;
}

cl_vector *
cl_matrix_transpose_dot_vector (const cl_matrix *a, const cl_vector *v)
{
	cl_vector *r;
	if (cl_matrix_is_empty (a)) cl_error ("cl_matrix_transpose_dot_vector", "matrix is empty");
	if (cl_vector_is_empty (v)) cl_error ("cl_matrix_transpose_dot_vector", "vector is empty");
	if (a->size1 != v->size) cl_error ("cl_matrix_transpose_dot_vector", "vector and matrix size not match");
	r = cl_vector_alloc (a->size2);
	if (cl_vector_is_empty (r)) cl_error ("cl_matrix_transpose_dot_vector", "vector not allocated");
	cblas_dgemv (CblasColMajor, CblasTrans, a->size1, a->size2, 1.0, a->data, a->lda, v->data, v->stride, 0.0, r->data, r->stride);
	return r;
}

