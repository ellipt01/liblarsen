/*
 * estimator.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdio.h>
#include <larsen.h>

#include "larsen_private.h"

/*
 *  Evaluate solution for L-1 (l->lambda2 = 0) or L-1 and L-n competitive (l->lambda2 != 0)
 *  penalized regression.
 *
 *  <input>
 *  larsen	*l		:	structure larsen.
 *  int	maxiter:	tolerance of the number of iterations for each regression steps
 *
 *  This function estimates regression coefficients of the system
 *
 *  beta = argmin | [l->y; 0] - l->scale * [l->x; sqrt(l->lambda2) * E] * beta |^2
 *  subject to | l->beta | <= l->lambda1
 *
 *  The beta that firstly satisfy lambda1 < sum |beta| is stored in l->beta.
 *  Its value (navie and elastic net solution) can be obtained by a function larsen_get_beta().
 */
bool
larsen_estimater (larsen *l, int maxiter)
{
	int		iter = 0;
	double	lambda1 = larsen_get_lambda1 (l, true);

	/* loop of regression */
	while (l->nrm1 <= lambda1 && !l->stop_loop) {
		if (!larsen_regression_step (l)) return false;
		if (++iter > maxiter) {
			fprintf (stderr, "number of iterations reaches max tolerance.\nregression stopped.\n");
			return false;
		}
	};

	/* when reached OLS but specified lambda1 is greater than |beta_ols|,
	 * stop regression */
	if (l->stop_loop && l->nrm1 < lambda1) return false;

	return true;
}
