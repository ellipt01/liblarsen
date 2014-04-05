/*
 * example.h
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#ifndef EXAMPLE_H_
#define EXAMPLE_H_

typedef struct {
	cl_matrix	*x;
	cl_vector	*y;
} larsen_data;

/* example.c */
larsen_data	*larsen_data_alloc (const size_t size1, const size_t size2);
void			larsen_data_free (larsen_data *data);
larsen_data	*read_data (char *fn, int skip_header);
void			output_solutionpath (int iter, larsen *l);
void			fprintf_beta (FILE *stream, int iter, larsen *l);

/* example_lasso.c */
void			example_elasticnet (cl_matrix *x0, cl_vector *y0, double start, double dt, double stop, double lambda2, double gamma, int maxsteps);

#endif /* EXAMPLE_H_ */
