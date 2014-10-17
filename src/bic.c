/*
 * bic.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdlib.h>
#include <math.h>
#include <larsen.h>

#include "private.h"

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
	double		*r = (double *) malloc (l->mu->nz * sizeof (double));
	double		*beta = larsen_copy_beta (l, true);	// scale * beta
	double		*mu = larsen_copy_mu (l, true);		// scale^2 * mu
	dcopy_ (&l->lreg->y->nz, l->lreg->y->data, &ione, r, &ione);
	daxpy_ (&l->mu->nz, &dmone, mu, &ione, r, &ione);	// r = - mu + r
	rss = pow (dnrm2_ (&l->mu->nz, r, &ione), 2.);
	if (!l->lreg->is_regtype_lasso) rss += l->lreg->lambda2 * pow (dnrm2_ (&l->beta->nz, beta, &ione), 2.);
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
	double	m = (double) l->lreg->x->m;
	double	n = (double) l->lreg->x->n;
	if (!l->lreg->is_regtype_lasso) m += n;

	return log (rss) + df * (log (m) + 2. * gamma * log (n)) / m;
}
