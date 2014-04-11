/*
 * bic.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <larsen.h>

/*
 *   Bayesian Information Criterion for L2 reguralized
 *   linear regression model b = Z * beta
 *   where b = [l->y ; 0], Z = [l->x ; sqrt(l->lambda2) * E]
 */

/* residual sum of squares | b - Z * beta |^2 */
static double
calc_rss (larsen *l)
{
	double		rss;
	double		*r = (double *) malloc (l->n * sizeof (double));
	double		*beta = larsen_get_beta (l);	// scale * beta
	double		*mu = larsen_get_mu (l);			// scale^2 * mu
	cblas_dcopy (l->n, l->y, 1, r, 1);
	cblas_daxpy (l->n, -1., mu, 1, r, 1);	// r = - mu + r
	rss = pow (cblas_dnrm2 (l->n, r, 1), 2.);
	if (l->is_elnet) rss += l->lambda2 * pow (cblas_dnrm2 (l->p, beta, 1), 2.);
	free (r);
	free (beta);
	free (mu);
	return rss;
}

/* degree of freedom
 * it is equal to #{j ; j \in A i.e., beta_j != 0} (Efron et al., 2004) */
static double
calc_degree_of_freedom (larsen *l)
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
larsen_eval_bic (larsen *l, double gamma)
{
	double	rss = calc_rss (l);
	double	df = calc_degree_of_freedom (l);
	double	n = (double) l->n;
	double	p = (double) l->p;
	if (l->is_elnet) n += p;

	return log (rss) + df * (log (n) + 2. * gamma * log (p)) / n;
}
