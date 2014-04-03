/*
 * cl_matrix.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include <cblas.h>

#include "cl_matrix.h"

void
cl_error (const char * function_name, const char *error_msg)
{
	fprintf (stderr, "ERROR: %s: %s\n", function_name, error_msg);
	exit (1);
}

int
get_index_of_matrix (const cl_matrix *a, int i, int j)
{
	return i + j * a->lda;
}

int
get_index_of_vector (const cl_vector *v, int i)
{
	return i * v->stride;
}

int
get_index_of_vector_int (const cl_vector_int *v, int i)
{
	return i * v->stride;
}

cl_matrix *
cl_matrix_alloc (const size_t size1, const size_t size2)
{
	cl_matrix	*a = (cl_matrix *) malloc (sizeof (cl_matrix));

	a->data = (double *) malloc (size1 * size2 * sizeof (double));
	if (!a->data) return NULL;

	a->size1 = size1;
	a->size2 = size2;
	a->lda = size1;
	a->tsize = a->lda * size2;
	a->owner = true;
	a->allocated = true;

	return a;
}

cl_vector *
cl_vector_alloc (const size_t size)
{
	cl_vector	*v = (cl_vector *) malloc (sizeof (cl_vector));

	v->data = (double *) malloc (size * sizeof (double));
	if (!v->data) return NULL;

	v->size = size;
	v->stride = 1;
	v->tsize = v->stride * v->size;
	v->owner = true;
	v->allocated = true;

	return v;
}

cl_vector_int *
cl_vector_int_alloc (const size_t size)
{
	cl_vector_int	*v = (cl_vector_int *) malloc (sizeof (cl_vector_int));

	v->data = (int *) malloc (size * sizeof (int));
	if (!v->data) return NULL;

	v->size = size;
	v->stride = 1;
	v->tsize = v->stride * v->size;
	v->owner = true;
	v->allocated = true;

	return v;
}

cl_matrix *
cl_matrix_view_array (const size_t size1, const size_t size2, double *data)
{
	cl_matrix *a = (cl_matrix *) malloc (sizeof (cl_matrix));
	a->size1 = size1;
	a->size2 = size2;
	a->lda = size1;
	a->tsize = size1 * size2;
	a->data = data;
	a->owner = false;
	a->allocated = true;
	return a;
}

bool
cl_matrix_is_empty (const cl_matrix *a)
{
	if (!a) return true;
	if (a->size1 <= 0 || a->size2 <= 0) return true;
	if (!a->allocated) return true;
	if (!a->data) return true;
	return false;
}

bool
cl_vector_is_empty (const cl_vector *v)
{
	if (!v) return true;
	if (v->size <= 0) return true;
	if (!v->allocated) return true;
	if (!v->data) return true;
	return false;
}

bool
cl_vector_int_is_empty (const cl_vector_int *v)
{
	if (!v) return true;
	if (v->size <= 0) return true;
	if (!v->allocated) return true;
	if (!v->data) return true;
	return false;
}

bool
cl_matrix_is_square (const cl_matrix *a)
{
	return (a->size1 == a->size2);
}

void
cl_matrix_free (cl_matrix *a)
{
	if (a && a->owner) {
		if (a->data) free (a->data);
		free (a);
	}
	return;
}

void
cl_vector_free (cl_vector *v)
{
	if (v && v->owner) {
		if (v->data) free (v->data);
		free (v);
	}
	return;
}

void
cl_vector_int_free (cl_vector_int *v)
{
	if (v && v->owner) {
		if (v->data) free (v->data);
		free (v);
	}
	return;
}

void
cl_matrix_set (cl_matrix *a, const int i, const int j, double val)
{
	int		index = get_index_of_matrix (a, i, j);
	if (index < 0 || a->tsize <= index) cl_error ("cl_matrix_set", "index out of range.");
	a->data[index] = val;
	return;
}

double
cl_matrix_get (const cl_matrix *a, const int i, const int j)
{
	int		index = get_index_of_matrix (a, i, j);
	if (index < 0 || a->tsize <= index) cl_error ("cl_matrix_get", "index out of range.");
	return a->data[index];
}

void
cl_vector_set (cl_vector *v, const int i, double val)
{
	int		index = get_index_of_vector (v, i);
	if (index < 0 || v->size * v->stride <= index) cl_error ("cl_vector_set", "index out of range.");
	v->data[index] = val;
	return;
}

double
cl_vector_get (const cl_vector *v, const int i)
{
	int		index = get_index_of_vector (v, i);
	if (index < 0 || v->tsize <= index) cl_error ("cl_vector_get", "index out of range.");
	return v->data[index];
}

void
cl_vector_int_set (cl_vector_int *v, const int i, int val)
{
	int		index = get_index_of_vector_int (v, i);
	if (index < 0 || v->size * v->stride <= index) cl_error ("cl_vector_set", "index out of range.");
	v->data[index] = val;
	return;
}

int
cl_vector_int_get (const cl_vector_int *v, const int i)
{
	int		index = get_index_of_vector_int (v, i);
	if (index < 0 || v->tsize <= index) cl_error ("cl_vector_get", "index out of range.");
	return v->data[index];
}

static cl_vector *
cl_vector_view_array (const size_t size, double *data)
{
	cl_vector *v = (cl_vector *) malloc (sizeof (cl_vector));
	v->size = size;
	v->stride = 1;
	v->data = data;
	v->owner = false;
	v->allocated = true;
	return v;
}

static cl_vector *
cl_vector_view_array_with_stride (const size_t size, const size_t stride, double *data)
{
	cl_vector	*v = cl_vector_view_array (size, data);
	v->stride = stride;
	return v;
}

cl_vector *
cl_matrix_column (const cl_matrix *a, const size_t k)
{
	cl_vector	*col;
	if (cl_matrix_is_empty (a)) cl_error ("cl_matrix_column", "input matrix is empty");
	if (k < 0 || a->size2 <= k) cl_error ("cl_matrix_column", "k must be 0 <= k < a->size2");
	col = cl_vector_view_array_with_stride (a->size1, 1, a->data + k * a->lda);
	return col;
}

void
cl_matrix_set_col (cl_matrix *a, const size_t index, const cl_vector *v)
{
	size_t incv, incm;

	if (cl_vector_is_empty (v)) cl_error ("cl_matrix_set_col", "vector is empty");
	if (cl_matrix_is_empty (a)) cl_error ("cl_matrix_get_col", "matrix is empty");
	if (v->size != a->size1) cl_error ("cl_matrix_set_col", "vector and matrix size invalid");
	if (index >= a->size2) cl_error ("cl_matrix_set_col", "specified index out of range");
	incv = v->stride;
	incm = 1;
	cblas_dcopy (v->size, v->data, incv, a->data + a->lda * index, incm);
	return;
}

void
cl_matrix_add_row (cl_matrix *a)
{
	size_t	lda;

	if (cl_matrix_is_empty (a)) cl_error ("cl_matrix_add_row", "matrix *a is empty");

	lda = a->lda;

	a->size1++;
	a->lda++;
	a->tsize = a->lda * a->size2;
	a->data = (double *) realloc (a->data, a->tsize * sizeof (double));

	{
		int			j;
		cl_vector	*col = cl_vector_alloc (a->size1);

		for (j = a->size2 - 1; 0 < j; j--) {
			cblas_dcopy (col->size, a->data + j * lda, 1, col->data, 1);
			cblas_dcopy (col->size, col->data, 1, a->data + j * a->lda, 1);
		}
		cl_vector_free (col);
		for (j = 0; j < a->size2; j++) cl_matrix_set (a, a->size1 - 1, j, 0.);
	}
	return;
}

void
cl_matrix_add_col (cl_matrix *a)
{
	size_t	lda;

	if (cl_matrix_is_empty (a)) cl_error ("cl_matrix_add_col", "matrix *a is empty");

	lda = a->lda;

	a->size2++;
	a->tsize = a->lda * a->size2;
	a->data = (double *) realloc (a->data, a->tsize * sizeof (double));

	return;
}

void
cl_matrix_remove_row (cl_matrix *a)
{
	size_t	lda;

	if (cl_matrix_is_empty (a)) cl_error ("cl_matrix_remove_row", "matrix *a is empty");

	lda = a->lda;

	a->size1--;
	a->lda--;
	a->tsize = a->lda * a->size2;

	if (a->size1 > 0) {
		int			j;
		cl_vector	*col = cl_vector_alloc (a->size1);
		for (j = 1; j < a->size2; j++) {
			cblas_dcopy (col->size, a->data + j * lda, 1, col->data, 1);
			cblas_dcopy (col->size, col->data, 1, a->data + j * a->lda, 1);
		}
		cl_vector_free (col);
	}
	if (a->data) a->data = (double *) realloc (a->data, a->tsize * sizeof (double));

	return;
}

void
cl_matrix_remove_col (cl_matrix *a)
{
	size_t	lda;

	if (cl_matrix_is_empty (a)) cl_error ("cl_matrix_remove_col", "matrix *a is empty");

	lda = a->lda;

	a->size2--;
	a->tsize = a->lda * a->size2;
	if (a->data) a->data = (double *) realloc (a->data, a->tsize * sizeof (double));

	return;
}

void
cl_matrix_memcpy (cl_matrix *dest, const cl_matrix *src)
{
	int		j;
	if (cl_matrix_is_empty (src)) cl_error ("cl_matrix_memcpy", "first matrix is empty");
	if (cl_matrix_is_empty (dest)) cl_error ("cl_matrix_memcpy", "second matrix is empty");
	if (dest->size1 != src->size1 || dest->size2 != src->size2) cl_error ("cl_matrix_memcpy", "matrix size not match");
	for (j = 0; j < src->size2; j++) {
		cl_vector	*col = cl_matrix_column ((cl_matrix *) src, j);
		cl_matrix_set_col (dest, j, col);
		cl_vector_free (col);
	}
	return;
}

void
cl_vector_memcpy (cl_vector *dest, const cl_vector *src)
{
	if (cl_vector_is_empty (src)) cl_error ("cl_vector_memcpy", "vector is empty");
	if (cl_vector_is_empty (dest)) cl_error ("cl_vector_memcpy", "vector not allocated");
	if (dest->size != src->size) cl_error ("cl_vector_memcpy", "vector size not match");
  	cblas_dcopy (src->size, src->data, src->stride, dest->data, dest->stride);
	return;
}

void
cl_matrix_set_zero (cl_matrix *a)
{
	int		i;
	if (cl_matrix_is_empty (a)) cl_error ("cl_matrix_set_zero", "matrix not allocated");
	for (i = 0; i < a->tsize; i++) a->data[i] = 0.0;
	return;
}

void
cl_vector_set_zero (cl_vector *v)
{
	int		i;
	if (cl_vector_is_empty (v)) cl_error ("cl_vector_set_zero", "vector not allocated");
	for (i = 0; i < v->size; i++) v->data[get_index_of_vector (v, i)] = 0.0;
	return;
}

void
cl_matrix_fprintf (FILE *stream, const cl_matrix *a, const char *format)
{
	int	i, j, k;
	if (cl_matrix_is_empty (a)) cl_error ("cl_matrix_fprintf2", "matrix is empty");
	k = 0;
	for (i = 0; i < a->size1; i++) {
		for (j = 0; j < a->size2; j++) {
			fprintf (stream, format, cl_matrix_get (a, i, j));
			if (j < a->size2 - 1) fprintf (stream, "%s", " ");
			k++;
		}
		fprintf (stream, "\n");
	}
	return;
}

void
cl_vector_fprintf (FILE *stream, const cl_vector *v, const char *format)
{
	int i;
	if (cl_vector_is_empty (v)) cl_error ("cl_vector_fprintf", "vector is empty");
	for (i = 0; i < v->size; i++) {
		fprintf (stream, format, v->data[get_index_of_vector (v, i)]);
		fprintf (stream, "\n");
	}
	return;
}
