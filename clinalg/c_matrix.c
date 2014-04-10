/*
 * c_matrix.c
 *
 *  Created on: 2014/04/03
 *      Author: utsugi
 */

#include "c_matrix.h"

/* blas */
extern void	dcopy_ (long *n, double *x, long *incx, double *y, long *incy);

void
c_error (const char * function_name, const char *error_msg)
{
	fprintf (stderr, "ERROR: %s: %s\n", function_name, error_msg);
	exit (1);
}

static c_matrix *
_allocate_c_matrix (void)
{
	c_matrix	*a = (c_matrix *) malloc (sizeof (c_matrix));
	if (!a) return NULL;
	a->size1 = 0;
	a->size2 = 0;
	a->lda = 0;
	a->tsize = 0;
	a->data = NULL;
	a->owner = false;
	return a;
}

static c_vector_int *
_allocate_c_vector_int (void)
{
	c_vector_int	*v = (c_vector_int *) malloc (sizeof (c_vector_int));
	if (!v) return NULL;
	v->size = 0;
	v->stride = 0;
	v->tsize = 0;
	v->data = NULL;
	v->owner = false;
	return v;
}

static c_vector *
_allocate_c_vector (void)
{
	c_vector	*v = (c_vector *) malloc (sizeof (c_vector));
	if (!v) return NULL;
	v->size = 0;
	v->stride = 0;
	v->tsize = 0;
	v->data = NULL;
	v->owner = false;
	return v;
}

c_matrix *
c_matrix_alloc (const size_t size1, const size_t size2)
{
	c_matrix	*a;

	if (size1 <= 0 || size2 <= 0) c_error ("c_matrix_alloc", "invalid size.");

	a = _allocate_c_matrix ();
	if (!a) {
		fprintf (stderr, "WARNING: c_matrix_alloc: failed to allocate memory.\n");
		return NULL;
	}
	a->data = (double *) malloc (size1 * size2 * sizeof (double));
	if (!a->data) fprintf (stderr, "WARNING: c_matrix_alloc: failed to allocate memory.\n");
	else {
		a->size1 = size1;
		a->size2 = size2;
		a->lda = size1;
		a->tsize = a->lda * a->size2;
		a->owner = true;
	}
	return a;
}

c_vector *
c_vector_alloc (const size_t size)
{
	c_vector	*v;

	if (size <= 0) c_error ("c_vector_alloc", "invalid size.");

	v = _allocate_c_vector ();
	if (!v) {
		fprintf (stderr, "WARNING: c_vector_alloc: failed to allocate memory.\n");
		return NULL;
	}
	v->data = (double *) malloc (size * sizeof (double));
	if (!v->data) fprintf (stderr, "WARNING: c_vector_alloc: failed to allocate memory.\n");
	else {
		v->size = size;
		v->stride = 1;
		v->tsize = v->stride * v->size;
		v->owner = true;
	}
	return v;
}

c_vector_int *
c_vector_int_alloc (const size_t size)
{
	c_vector_int	*v;

	if (size <= 0) c_error ("c_vector_int_alloc", "invalid size.");

	v = _allocate_c_vector_int ();
	if (!v) {
		fprintf (stderr, "WARNING: c_vector_int_alloc: failed to allocate memory.\n");
		return NULL;
	}
	v->data = (int *) malloc (size * sizeof (int));
	if (!v->data) fprintf (stderr, "WARNING: c_vector_int_alloc: failed to allocate memory.\n");
	else {
		v->size = size;
		v->stride = 1;
		v->tsize = v->stride * v->size;
		v->owner = true;
	}
	return v;
}

c_matrix *
c_matrix_view_array (const size_t size1, const size_t size2, const size_t lda, double *data)
{
	c_matrix *a;

	if (!data) c_error ("c_matrix_view_array", "array is empty.");
	if (size1 <= 0 || size2 <= 0 || lda <= 0) c_error ("c_matrix_view_array", "invalid size.");

	a = _allocate_c_matrix ();
	a->size1 = size1;
	a->size2 = size2;
	a->lda = lda;
	a->tsize = a->lda * a->size2;
	a->data = data;
	a->owner = false;
	return a;
}

c_vector *
c_vector_view_array (const size_t size, const size_t stride, double *data)
{
	c_vector *v;

	if (!data) c_error ("c_vector_view_array", "array is empty.");
	if (size <= 0 || stride <= 0) c_error ("c_vector_view_array", "invalid size.");

	v = _allocate_c_vector ();
	v->size = size;
	v->stride = stride;
	v->tsize = v->stride * v->size;
	v->data = data;
	v->owner = false;
	return v;
}

bool
c_matrix_is_empty (const c_matrix *a)
{
	if (!a) return true;
	if (!a->data) return true;
	return false;
}

bool
c_vector_is_empty (const c_vector *v)
{
	if (!v) return true;
	if (!v->data) return true;
	return false;
}

bool
c_vector_int_is_empty (const c_vector_int *v)
{
	if (!v) return true;
	if (!v->data) return true;
	return false;
}

bool
c_matrix_is_square (const c_matrix *a)
{
	return (a->size1 == a->size2);
}

void
c_matrix_free (c_matrix *a)
{
	if (a) {
		if (a->data && a->owner) free (a->data);
		free (a);
	}
	return;
}

void
c_vector_free (c_vector *v)
{
	if (v) {
		if (v->data && v->owner) free (v->data);
		free (v);
	}
	return;
}

void
c_vector_int_free (c_vector_int *v)
{
	if (v) {
		if (v->data && v->owner) free (v->data);
		free (v);
	}
	return;
}

void
c_matrix_set (c_matrix *a, const int i, const int j, double val)
{
	int		index;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_set", "matrix is empty.");
	index = GET_INDEX_OF_MATRIX (a, i, j);
	if (index < 0 || a->tsize <= index) c_error ("c_matrix_set", "index out of range.");
	a->data[index] = val;
	return;
}

double
c_matrix_get (const c_matrix *a, const int i, const int j)
{
	int		index;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_get", "matrix is empty.");
	index = GET_INDEX_OF_MATRIX (a, i, j);
	if (index < 0 || a->tsize <= index) c_error ("c_matrix_get", "index out of range.");
	return a->data[index];
}

void
c_vector_set (c_vector *v, const int i, double val)
{
	int		index;
	if (c_vector_is_empty (v)) c_error ("c_vector_set", "vector is empty.");
	index = GET_INDEX_OF_VECTOR (v, i);
	if (index < 0 || v->size * v->stride <= index) c_error ("c_vector_set", "index out of range.");
	v->data[index] = val;
	return;
}

double
c_vector_get (const c_vector *v, const int i)
{
	int		index;
	if (c_vector_is_empty (v)) c_error ("c_vector_get", "vector is empty.");
	index = GET_INDEX_OF_VECTOR (v, i);
	if (index < 0 || v->tsize <= index) c_error ("c_vector_get", "index out of range.");
	return v->data[index];
}

void
c_vector_int_set (c_vector_int *v, const int i, int val)
{
	int		index;
	if (c_vector_int_is_empty (v)) c_error ("c_vector_int_set", "vector is empty.");
	index = GET_INDEX_OF_VECTOR (v, i);
	if (index < 0 || v->size * v->stride <= index) c_error ("c_vector_int_set", "index out of range.");
	v->data[index] = val;
	return;
}

int
c_vector_int_get (const c_vector_int *v, const int i)
{
	int		index;
	if (c_vector_int_is_empty (v)) c_error ("c_vector_int_get", "vector is empty.");
	index = GET_INDEX_OF_VECTOR (v, i);
	if (index < 0 || v->tsize <= index) c_error ("c_vector_int_get", "index out of range.");
	return v->data[index];
}

void
c_matrix_get_col (c_vector *y, const c_matrix *a, const size_t index)
{
	long	n;
	long	incx;
	long	incy;
	if (c_vector_is_empty (y)) c_error ("c_matrix_get_col", "vector is empty.");
	if (c_matrix_is_empty (a)) c_error ("c_matrix_get_col", "matrix is empty.");
	if (y->size != a->size1) c_error ("c_matrix_get_col", "vector and matrix size does not match.");
	if (index < 0 || index >= a->size2) c_error ("c_matrix_get_col", "index must be in [0, a->size2).");
	n = (long) a->size1;
	incx = 1;
	incy = (long) y->stride;
	dcopy_ (&n, a->data + a->lda * index, &incx, y->data, &incy);
	return;
}

void
c_matrix_set_col (c_matrix *a, const size_t index, const c_vector *x)
{
	long	n;
	long	incx;
	long	incy;
	if (c_vector_is_empty (x)) c_error ("c_matrix_set_col", "vector is empty.");
	if (c_matrix_is_empty (a)) c_error ("c_matrix_set_col", "matrix is empty.");
	if (x->size != a->size1) c_error ("c_matrix_set_col", "vector and matrix size does not match.");
	if (index < 0 || index >= a->size2) c_error ("c_matrix_set_col", "index must be in [0, a->size2).");
	n = (long) x->size;
	incx = (long) x->stride;
	incy = 1;
	dcopy_ (&n, x->data, &incx, a->data + a->lda * index, &incy);
	return;
}

void
c_matrix_add_row (c_matrix *a)
{
	long	n;
	long	incx = 1;
	long	incy = 1;
	size_t	lda;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_add_row", "matrix is empty.");

	lda = a->lda;
	a->size1++;
	a->lda++;
	a->tsize = a->lda * a->size2;
	a->data = (double *) realloc (a->data, a->tsize * sizeof (double));

	{
		int			j;
		c_vector	*col = c_vector_alloc (a->size1);

		n = (long) a->size1;
		for (j = a->size2 - 1; 0 < j; j--) {
			dcopy_ (&n, a->data + j * lda, &incx, col->data, &incy);
			dcopy_ (&n, col->data, &incy, a->data + j * a->lda, &incx);
		}
		c_vector_free (col);
		for (j = 0; j < a->size2; j++) c_matrix_set (a, a->size1 - 1, j, 0.);
	}
	return;
}

void
c_matrix_add_col (c_matrix *a)
{
	size_t	lda;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_add_col", "matrix is empty.");

	lda = a->lda;
	a->size2++;
	a->tsize = a->lda * a->size2;
	a->data = (double *) realloc (a->data, a->tsize * sizeof (double));

	return;
}

void
c_matrix_remove_row (c_matrix *a)
{
	long	n;
	long	incx = 1;
	long	incy = 1;
	size_t	lda;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_remove_row", "matrix is empty.");

	lda = a->lda;
	a->size1--;
	a->lda--;
	a->tsize = a->lda * a->size2;
	n = (long) a->size2;
	if (a->size1 > 0) {
		int			j;
		c_vector	*col = c_vector_alloc (a->size1);
		for (j = 1; j < a->size2; j++) {
			dcopy_ (&n, a->data + j * lda, &incx, col->data, &incy);
			dcopy_ (&n, col->data, &incx, a->data + j * a->lda, &incx);
		}
		c_vector_free (col);
	}
	if (a->data) a->data = (double *) realloc (a->data, a->tsize * sizeof (double));

	return;
}

void
c_matrix_remove_col (c_matrix *a)
{
	size_t	lda;

	if (c_matrix_is_empty (a)) c_error ("c_matrix_remove_col", "matrix is empty.");

	lda = a->lda;
	a->size2--;
	a->tsize = a->lda * a->size2;
	if (a->data) a->data = (double *) realloc (a->data, a->tsize * sizeof (double));

	return;
}

void
c_matrix_memcpy (c_matrix *dest, const c_matrix *src)
{
	int			j;
	c_vector	*col;
	if (c_matrix_is_empty (src)) c_error ("c_matrix_memcpy", "first matrix is empty.");
	if (c_matrix_is_empty (dest)) c_error ("c_matrix_memcpy", "second matrix is empty.");
	if (dest->size1 != src->size1 || dest->size2 != src->size2) c_error ("c_matrix_memcpy", "matrix size does not match.");
	col = c_vector_alloc (src->size1);
	for (j = 0; j < src->size2; j++) {
		c_matrix_get_col (col, src, j);
		c_matrix_set_col (dest, j, col);
	}
	c_vector_free (col);
	return;
}

void
c_vector_memcpy (c_vector *dest, const c_vector *src)
{
	long	n;
	long	incx;
	long	incy;
	if (c_vector_is_empty (src)) c_error ("c_vector_memcpy", "first vector is empty.");
	if (c_vector_is_empty (dest)) c_error ("c_vector_memcpy", "second vector is empty.");
	if (dest->size != src->size) c_error ("c_vector_memcpy", "vector size does not match.");
	n = (long) src->size;
	incx = (long) src->stride;
	incy = (long) dest->stride;
	dcopy_ (&n, src->data, &incx, dest->data, &incy);
	return;
}

void
c_matrix_set_zero (c_matrix *a)
{
	int		i;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_set_zero", "matrix is empty.");
	for (i = 0; i < a->tsize; i++) a->data[i] = 0.0;
	return;
}

void
c_vector_set_zero (c_vector *v)
{
	int		i;
	if (c_vector_is_empty (v)) c_error ("c_vector_set_zero", "vector is empty.");
	for (i = 0; i < v->size; i++) v->data[GET_INDEX_OF_VECTOR (v, i)] = 0.0;
	return;
}

void
c_matrix_fprintf (FILE *stream, const c_matrix *a, const char *format)
{
	int	i, j, k;
	if (c_matrix_is_empty (a)) c_error ("c_matrix_fprintf", "matrix is empty.");
	k = 0;
	for (i = 0; i < a->size1; i++) {
		for (j = 0; j < a->size2; j++) {
			fprintf (stream, format, c_matrix_get (a, i, j));
			if (j < a->size2 - 1) fprintf (stream, "%s", " ");
			k++;
		}
		fprintf (stream, "\n");
	}
	return;
}

void
c_vector_fprintf (FILE *stream, const c_vector *v, const char *format)
{
	int i;
	if (c_vector_is_empty (v)) c_error ("c_vector_fprintf", "vector is empty.");
	for (i = 0; i < v->size; i++) {
		fprintf (stream, format, v->data[GET_INDEX_OF_VECTOR (v, i)]);
		fprintf (stream, "\n");
	}
	return;
}
