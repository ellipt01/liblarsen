/*
 * bic.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdlib.h>
#include <math.h>
#include <larsen.h>

#include "larsen_private.h"

/*
 *   Bayesian Information Criterion for L2 reguralized
 *   linear regression model b = Z * beta
 *   where b = [l->y ; 0], Z = [l->x ; sqrt(l->lambda2) * E]
 */

/* residual sum of squares | b - Z * beta |^2 */
static double
calc_rss (const larsen *l)
{
	size_t			n = linsys_get_n (l->lsys);
	size_t			p = linsys_get_p (l->lsys);
	const double	*y = linsys_get_y (l->lsys);
	double			lambda2 = linsys_get_lambda2 (l->lsys);

	double			rss;
	double			*r = (double *) malloc (n * sizeof (double));
	double			*beta;
	double			*mu;

	dcopy_ (LINSYS_CINTP (n), y, &ione, r, &ione);
	mu = larsen_copy_mu (l, true);	// = scale^2 * mu
	daxpy_ (LINSYS_CINTP (n), &dmone, mu, &ione, r, &ione);	// r = - mu + r
	free (mu);

	rss = pow (dnrm2_ (LINSYS_CINTP (n), r, &ione), 2.);
	free (r);

	beta = larsen_copy_beta (l, true);	// = scale * beta
	if (!linsys_is_regtype_lasso (l->lsys)) rss += lambda2 * pow (dnrm2_ (LINSYS_CINTP (p), beta, &ione), 2.);
	free (beta);

	return rss;
}

/* degree of freedom
 * it is equal to #{j ; j \in A i.e., beta_j != 0} (Efron et al., 2004) */
static double
calc_degree_of_freedom (const larsen *l)
{
	return (double) l->sizeA;
}

/* Extended Bayesian Information Criterion (Chen and Chen, 2008)
 * EBIC = n log(rss) + df * log(n) + 2 * gamma * df * log(p)
 * gamma	: tuning parameter for EBIC
 * rss		: residual sum of squares |b - Z * beta|^2
 * df		: degree of freedom of the system
 * 	n		: number of data (= l->y->size = l->x->size1)
 * 	p		: number of variables (= l->x->size2)
 *
 * 	if gamma = 0, eBIC is identical with the classical BIC
*/
double
larsen_eval_bic (const larsen *l, double gamma)
{
	double	rss = calc_rss (l);
	double	df = calc_degree_of_freedom (l);
	double	n = (double) linsys_get_n (l->lsys);
	double	p = (double) linsys_get_p (l->lsys);
	if (!linsys_is_regtype_lasso (l->lsys)) n += p;

	return log (rss) + df * (log (n) + 2. * gamma * log (p)) / n;
}
