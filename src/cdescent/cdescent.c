/*
 * cdescent.c
 *
 *  Created on: 2014/05/27
 *      Author: utsugi
 */

#include <stdlib.h>
#include <math.h>
#include <cdescent.h>

#include "linreg_private.h"

static double
soft_threshold (const double z, const double gamma)
{
	double	val = 0.;
	if (gamma < fabs (z)) val = (z > 0.) ? z - gamma : z + gamma;
	return val;
}

bool
cdescent_cyclick_step (cdescent *cd)
{
	int				j;
	size_t			n = cd->lreg->n;
	size_t			p = cd->lreg->p;
	const double	*x = cd->lreg->x;
	double			lambda1 = cd->lambda1;
	double			scale2 = cd->lreg->scale2;

	double			nrm;
	double			*delta_beta = (double *) malloc (p * sizeof (double));
	bool			converged;

	// penalty
	size_t			pj = cd->lreg->pen->pj;
	const double	*jr = cd->lreg->pen->r;
	double			*jb = NULL;	// J * beta
	if (cdescent_is_regtype_userdef (cd)) {
		// J * beta
		jb = (double *) malloc (pj * sizeof (double));
		dgemv_ ("N", LINREG_CINTP (pj), LINREG_CINTP (p), &done, jr, LINREG_CINTP (pj), cd->beta, &ione, &dzero, jb, &ione);
	}

	cd->nrm1 = cd->nrm1_prev;
	dcopy_ (LINREG_CINTP (p), cd->beta, &ione, cd->beta_prev, &ione);

	for (j = 0; j < p; j++) {
		/*** z = X(:,j)' * y - X(:,j)' * X * beta + beta(j) ***/
		double			cj = cd->c[j];	// X' * y
		const double	*xj = x + LINREG_INDEX_OF_MATRIX (0, j, n);	// X(:,j)
		double			xjtm = ddot_ (LINREG_CINTP (n), xj, &ione, cd->mu, &ione);	// X(:,j)' * mu

		double			z = cj - xjtm + cd->beta[j];

		/* user defined penalty (not lasso and not ridge) */
		if (cdescent_is_regtype_userdef (cd)) {
			/*** z = z - lambda2 * J(:,j)' * J * beta
			 *         + lambda2 * J(:,j)' * J(:,j) * beta(j) ***/
			double			lambda2 = cd->lreg->lambda2;
			const double	*jrj = jr + LINREG_INDEX_OF_MATRIX (0, j, pj);	// J(:,j)
			// J(:,j)' * J(:,j)
			double			jtj = pow (dnrm2_ (LINREG_CINTP (pj), jrj, &ione), 2.);
			// z -= lambda2 * J(:,j)' * (J * beta)
			z -= lambda2 * ddot_ (LINREG_CINTP (pj), jrj, &ione, jb, &ione);
			// z += lambda2 * J(:,j)' * J(:,j) * beta(j)
			z += lambda2 * jtj * cd->beta[j];
			scale2 = 1. / (1. + jtj * lambda2);
		}

		cd->beta[j] = scale2 * soft_threshold (z, lambda1);

		delta_beta[j] = cd->beta[j] - cd->beta_prev[j];
		if (fabs (delta_beta[j]) > 0.) {
			// mu += X(:, j) * (beta[j] - beta_prev[j])
			daxpy_ (LINREG_CINTP (n), &delta_beta[j], xj, &ione, cd->mu, &ione);

			if (cdescent_is_regtype_userdef (cd)) {
				// jb += J(:,j) * (beta[j] - beta_prev[j])
				size_t			pj = cd->lreg->pen->pj;
				const double	*jrj = cd->lreg->pen->r + LINREG_INDEX_OF_MATRIX (0, j, pj);	// J(:,j)
				daxpy_ (LINREG_CINTP (pj), &delta_beta[j], jrj, &ione, jb, &ione);
			}
		}

	}
	if (jb) free (jb);

	cd->nrm1 = dasum_ (LINREG_CINTP (p), cd->beta, &ione);
	nrm = dnrm2_ (LINREG_CINTP (p), delta_beta, &ione);
	free (delta_beta);

	converged = (nrm <= cd->tolerance);

	return converged;
}

bool
cdescent_cyclick (cdescent *cd, const int maxiter)
{
	int		iter = 0;
	bool	converged = false;

	while (!converged) {

		converged = cdescent_cyclick_step (cd);

		if (++iter >= maxiter) break;
	}

	return converged;
}
