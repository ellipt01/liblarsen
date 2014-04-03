/*
 * elasticnet.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include "larsen.h"

#define _DUMMY_ -1
extern double	larsen_get_lambda1 (larsen *l);

bool
larsen_elasticnet (larsen *l, int maxiter)
{
	int		iter = 0;
	double	lambda1 = larsen_get_lambda1 (l);
	double	nrm1 = cl_vector_asum (l->beta);
	while (nrm1 <= lambda1 && larsen_loop_continue (l, _DUMMY_)) {
		if (!larsen_regression_step (l)) return false;
		nrm1 = cl_vector_asum (l->beta);
		if (++iter > maxiter) {
			fprintf (stderr, "number of iterations reaches max tolerance.\nregression stopped.\n");
			return false;
		}
	};
	return larsen_interpolate (l);
}
