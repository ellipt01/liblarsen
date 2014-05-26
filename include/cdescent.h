/*
 * cdescent.h
 *
 *  Created on: 2014/05/27
 *      Author: utsugi
 */

#ifndef CDESCENT_H_
#define CDESCENT_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <linreg.h>

typedef struct s_cdescent	cdescent;

struct s_cdescent {

	const linreg	*lreg;

	double			tolerance;
	double			lambda1;

	double			*c;

	double			nrm1;
	double			*beta;
	double			*mu;

	double			nrm1_prev;
	double			*beta_prev;

};

/* utils.c */
cdescent	*cdescent_alloc (const linreg *lreg, const double lambda1, const double tol);
void		cdescent_free (cdescent *cd);

void		cdescent_set_lambda1 (cdescent *cd, const double lambda1);
double		cdescent_get_lambda2 (const cdescent *cd, bool scaling);

double		*cdescent_copy_beta (const cdescent *cd, bool scaling);
double		*cdescent_copy_mu (const cdescent *cd, bool scaling);

bool		cdescent_is_regtype_lasso (const cdescent *cd);
bool		cdescent_is_regtype_ridge (const cdescent *cd);
bool		cdescent_is_regtype_userdef (const cdescent *cd);

/* cdescent.c */
bool		cdescent_cyclick_step (cdescent *cd);
bool		cdescent_cyclick (cdescent *cd, const int maxiter);

#ifdef __cplusplus
}
#endif

#endif /* CDESCENT_H_ */
