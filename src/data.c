/*
 * data.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <larsen.h>

/* centering vector i.e., y -> y - mean(y) */
double
larsen_centering_vector (cl_vector *y)
{
	double	mean = cl_vector_mean (y);
	cl_vector_add_constant (y, - mean);
	return mean;
}

/* normalizing vector i.e., y -> y / norm(y) */
double
larsen_normalizing_vector (cl_vector *y)
{
	double	nrm = cl_vector_nrm (y);
	cl_vector_scale (y, 1. / nrm);
	return nrm;
}

/* centering each column vector of matrix */
cl_vector *
larsen_centering_matrix (cl_matrix *x)
{
	int			j;
	cl_vector	*mean = cl_vector_alloc (x->size2);
	for (j = 0; j < x->size2; j++) {
		cl_vector	*xj = cl_matrix_column (x, j);
		double		meanj = cl_vector_mean (xj);
		cl_vector_set (mean, j, meanj);
		cl_vector_add_constant (xj, - meanj);
	}
	return mean;
}

/* normalizing each column vector of matrix */
cl_vector *
larsen_normalizing_matrix (cl_matrix *x)
{
	int			j;
	cl_vector	*nrm = cl_vector_alloc (x->size2);
	for (j = 0; j < x->size2; j++) {
		cl_vector	*xj = cl_matrix_column (x, j);
		double		nrmj = cl_vector_nrm (xj);
		cl_vector_set (nrm, j, nrmj);
		cl_vector_scale (xj, 1. / nrmj);
	}
	return nrm;
}
