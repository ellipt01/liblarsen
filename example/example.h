/*
 * example.h
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#ifndef EXAMPLE_H_
#define EXAMPLE_H_

/* example.c */
void	read_data (char *fn, int skip_header, size_t *n, size_t *p, double **y, double **x);
void	output_solutionpath (int iter, larsen *l);
void	fprintf_beta (FILE *stream, int iter, larsen *l);

/* example_estimator.c */
void	example_l1regression (const linreg *lreg, double start, double dt, double stop, double gamma, int maxiter);

#endif /* EXAMPLE_H_ */
