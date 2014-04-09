/*
 * data.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <larsen.h>

/* centering vector: y -> y - mean(y) */
double
larsen_centering_vector (c_vector *y)
{
	double	mean = c_vector_mean (y);
	c_vector_add_constant (y, - mean);
	return mean;
}

/* normalizing vector: y -> y / norm(y) */
double
larsen_normalizing_vector (c_vector *y)
{
	double	nrm = c_vector_nrm (y);
	c_vector_scale (y, 1. / nrm);
	return nrm;
}

/* centering each column of matrix:
 * x(:, j) ->  x(:, j) - mean(x(:, j)) */
c_vector *
larsen_centering_matrix (c_matrix *x)
{
	int			j;
	c_vector	*xj = c_vector_alloc (x->size1);
	c_vector	*mean = c_vector_alloc (x->size2);
	for (j = 0; j < x->size2; j++) {
		double		meanj;
		c_matrix_get_col (xj, x, j);
		meanj = c_vector_mean (xj);
		c_vector_set (mean, j, meanj);
		c_vector_add_constant (xj, - meanj);
		c_matrix_set_col (x, j, xj);
	}
	c_vector_free (xj);
	return mean;
}

/* normalizing each column of matrix:
 * x(:, j) -> x(:, j) / norm(x(:, j)) */
c_vector *
larsen_normalizing_matrix (c_matrix *x)
{
	int			j;
	c_vector	*xj = c_vector_alloc (x->size1);
	c_vector	*nrm = c_vector_alloc (x->size2);
	for (j = 0; j < x->size2; j++) {
		double		nrmj;
		c_matrix_get_col (xj, x, j);
		nrmj = c_vector_nrm (xj);
		c_vector_set (nrm, j, nrmj);
		c_vector_scale (xj, 1. / nrmj);
		c_matrix_set_col (x, j, xj);
	}
	c_vector_free (xj);
	return nrm;
}
