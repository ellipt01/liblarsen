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
	cl_vector	*r = cl_vector_alloc (l->y->size);
	cl_vector	*beta = larsen_get_beta (l);
	cl_vector	*mu = larsen_get_mu (l);
	cl_vector_memcpy (r, l->y);
	cl_vector_sub (r, mu);
	rss = pow (cl_vector_nrm (r), 2.);
	if (l->do_scaling) rss += l->lambda2 * pow (cl_vector_nrm (beta), 2.);
	cl_vector_free (r);
	cl_vector_free (beta);
	cl_vector_free (mu);
	return rss;
}

/* degree of freedom
 * it is equal to #{j ; j \in A i.e., beta_j != 0} (Efron et al., 2004) */
static double
calc_degree_of_freedom (larsen *l)
{
	return (double) l->A->size;
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
	double	n = (double) l->y->size;
	double	p = (double) l->x->size2;
	if (l->do_scaling) n += p;

	return log (rss) + df * (log (n) + 2. * gamma * log (p)) / n;
}
