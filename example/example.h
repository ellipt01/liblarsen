/*
 * example.h
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#ifndef EXAMPLE_H_
#define EXAMPLE_H_

/* example.c */
void			read_data (char *fn, int skip_header, cl_vector **y, cl_matrix **x);
void			output_solutionpath (int iter, larsen *l);
void			fprintf_beta (FILE *stream, int iter, larsen *l);

/* example_elasticnet.c */
void			example_elasticnet (cl_matrix *x0, cl_vector *y0, double start, double dt, double stop, double lambda2, double gamma, int maxsteps);

#endif /* EXAMPLE_H_ */
