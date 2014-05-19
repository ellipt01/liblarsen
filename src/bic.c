/*
 * bic.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdlib.h>
#include <math.h>
#include <larsen.h>

#include "linsys_private.h"

/*
 *   Bayesian Information Criterion for L2 reguralized
 *   linear regression model b = Z * beta
 *   where b = [l->y ; 0], Z = [l->x ; sqrt(l->lambda2) * E]
 */

/* residual sum of squares | b - Z * beta |^2 */
static double
calc_rss (const larsen *l)
{
	double		rss;
	double		*r = (double *) malloc (l->lsys->n * sizeof (double));
	double		*beta = larsen_copy_beta_elasticnet (l);	// scale * beta
	double		*mu = larsen_copy_mu_elasticnet (l);		// scale^2 * mu
	dcopy_ (LINSYS_CINTP (l->lsys->n), l->lsys->y, &ione, r, &ione);
	daxpy_ (LINSYS_CINTP (l->lsys->n), &dmone, mu, &ione, r, &ione);	// r = - mu + r
	rss = pow (dnrm2_ (LINSYS_CINTP (l->lsys->n), r, &ione), 2.);
	if (!l->is_lasso) rss += l->lambda2 * pow (dnrm2_ (LINSYS_CINTP (l->lsys->p), beta, &ione), 2.);
	free (r);
	free (beta);
	free (mu);
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
	double	n = (double) l->lsys->n;
	double	p = (double) l->lsys->p;
	if (!l->is_lasso) n += p;

	return log (rss) + df * (log (n) + 2. * gamma * log (p)) / n;
}
