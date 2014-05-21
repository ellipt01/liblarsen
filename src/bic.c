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
	if (!linsys_is_regtype_lasso (l->lsys)) {
		if (linsys_is_regtype_ridge (l->lsys))
			// rss += lambda2 * |beta|^2
			rss += lambda2 * pow (dnrm2_ (LINSYS_CINTP (p), beta, &ione), 2.);
		else {
			size_t			pj = linsys_get_pj (l->lsys);
			const double	*jr = linsys_get_penalty (l->lsys);
			double			*jb = (double *) malloc (pj * sizeof (double));
			// J * beta
			dgemv_ ("N", LINSYS_CINTP (pj), LINSYS_CINTP (p), &done, jr, LINSYS_CINTP (pj), beta, &ione, &dzero, jb, &ione);
			// rss += lambda2 * |J * beta|^2
			rss += lambda2 * pow (dnrm2_ (LINSYS_CINTP (p), jb, &ione), 2.);
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
	double	n = (double) linsys_get_n (l->lsys);
	double	p = (double) linsys_get_p (l->lsys);
	if (!linsys_is_regtype_lasso (l->lsys)) n += p;

	return log (rss) + df * (log (n) + 2. * gamma * log (p)) / n;
}
