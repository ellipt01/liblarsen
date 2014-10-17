/*
 * example.h
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#ifndef EXAMPLE_H_
#define EXAMPLE_H_

/* example.c */
void	output_solutionpath (int iter, larsen *l);
void	fprintf_beta (FILE *stream, int iter, larsen *l);

/* example_elasticnet.c */
void	example_elasticnet (const linregmodel *lreg, double start, double dt, double stop, double gamma, int maxiter);

#endif /* EXAMPLE_H_ */
