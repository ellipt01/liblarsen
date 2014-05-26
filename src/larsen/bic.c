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
 *   Bayesian Information Criterion for a linear regression model
 *   b = Z * beta
 *   where b = [y ; 0], Z = [X ; sqrt(l->lambda2) * J]
 */

/* residual sum of squares | b - Z * beta |^2 */
static double
calc_rss (const larsen *l)
{
	size_t			n = l->lreg->n;
	size_t			p = l->lreg->p;
	const double	*y = l->lreg->y;
	double			lambda2 = l->lreg->lambda2;

	double			rss;
	double			*r = (double *) malloc (n * sizeof (double));
	double			*beta;
	double			*mu;

	dcopy_ (LINREG_CINTP (n), y, &ione, r, &ione);
	mu = larsen_copy_mu (l, true);	// = scale^2 * mu
	daxpy_ (LINREG_CINTP (n), &dmone, mu, &ione, r, &ione);	// r = - mu + r
	free (mu);

	rss = pow (dnrm2_ (LINREG_CINTP (n), r, &ione), 2.);
	free (r);

	beta = larsen_copy_beta (l, true);	// = scale * beta
	if (!larsen_is_regtype_lasso (l)) {
		if (larsen_is_regtype_ridge (l))
			// rss += lambda2 * |beta|^2
			rss += lambda2 * pow (dnrm2_ (LINREG_CINTP (p), beta, &ione), 2.);
		else {
			size_t			pj = l->lreg->pen->pj;
			const double	*jr = l->lreg->pen->r;
			double			*jb = (double *) malloc (pj * sizeof (double));
			// J * beta
			dgemv_ ("N", LINREG_CINTP (pj), LINREG_CINTP (p), &done, jr, LINREG_CINTP (pj), beta, &ione, &dzero, jb, &ione);
			// rss += lambda2 * |J * beta|^2
			rss += lambda2 * pow (dnrm2_ (LINREG_CINTP (p), jb, &ione), 2.);
			free (jb);
		}
	}
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
	double	n = (double) l->lreg->n;
	double	p = (double) l->lreg->p;
	if (!larsen_is_regtype_lasso (l)) n += p;

	return log (rss) + df * (log (n) + 2. * gamma * log (p)) / n;
}
