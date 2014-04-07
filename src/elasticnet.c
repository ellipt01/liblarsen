/*
 * elasticnet.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <larsen.h>

#define _DUMMY_ -1
extern double	larsen_get_lambda1 (larsen *l);

/*
 *  Evaluate lasso (l->lambda2 = 0) or elastic net (l->lambda2 != 0) estimator.
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
 *  The optimal beta corresponding to a designed l->lambda1 is stored
 *  if l->interp == false, in l->beta, else in l->beta_interp.
 *  Its value (elastic net solution) can be obtained by a function larsen_get_beta (larsen *l).
 */
bool
larsen_elasticnet (larsen *l, int maxiter)
{
	int		iter = 0;
	double	lambda1 = larsen_get_lambda1 (l);
	double	nrm1 = cl_vector_asum (l->beta);
	while (nrm1 <= lambda1 && !l->stop_loop) {
		if (!larsen_regression_step (l)) return false;
		nrm1 = cl_vector_asum (l->beta);
		if (++iter > maxiter) {
			fprintf (stderr, "number of iterations reaches max tolerance.\nregression stopped.\n");
			return false;
		}
	};

	// if reached OLS but specified lambda1 is greater than |beta_ols|, stop regression
	if (l->stop_loop && nrm1 < lambda1) return false;

	larsen_interpolate (l);
	return true;
}
