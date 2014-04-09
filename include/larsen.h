/*
 ============================================================================
 Name        : larsen.h
 Author      : Mitsuru Utsugi
 Version     :
 ============================================================================
 */

#ifndef LARSEN_H_
#define LARSEN_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include <c_linalg.h>

typedef enum {
	ACTIVESET_ACTION_NONE	= -1,
	ACTIVESET_ACTION_ADD		=  0,	// add new variable to the active set
	ACTIVESET_ACTION_DROP	=  1,	// remove a variable from the active set
	NUM_ACTIVESET_ACTION		=  2
} ActiveSetAction;

typedef struct {
	ActiveSetAction	action;
	int					column;	// index of operand column of matrix X
	int					index;		// position of added / removed item in A
} activeset_oper;

typedef struct s_larsen	larsen;

struct s_larsen {

	/* if true, loop of lasso or elastic net regression is terminated */
	bool					stop_loop;

	/* threshold for L1 penalty */
	double					lambda1;
	/* threshold for L2 penalty */
	double					lambda2;

	/* true: elastic net, false: lasso */
	bool					is_elnet;
	/* scale = (is_elnet) ? sqrt(1 + lambda2) : 1 */
	double					scale;

	const c_matrix		*x;		// predictors
	const c_vector		*y;		// response

	/* correlation */
	double					sup_c;	// sup(abs(c))
	c_vector				*c;		// colleration: = X' * (y - mu)

	/* active set */
	activeset_oper		oper;
	c_vector_int			*A;		// active set
	c_vector_int			*Ac;	// implement of A

	/* equiangular */
	double					absA;
	c_vector				*u;		// equiangular vector
	c_vector				*w;

	/* step_size */
	double					stepsize;	// step size which beta will be progress.

	/* solution */
	c_vector				*beta;	// regression coefficients
	c_vector				*mu;	// estimated response

	/* previous beta and mu */
	c_vector				*beta_prev;	// previous beta
	c_vector				*mu_prev;		// previous mu

	/* interpolation */
	bool					interp;		// interpolation was done or not.
	double					stepsize_intr;
	c_vector				*beta_intr;	// interpolated beta
	c_vector				*mu_intr;		// interpolated mu

	/* cholesky factorization */
	c_matrix				*chol;	// = chol(XA' * XA), where XA = X(A).

};

/* util.c */
larsen		*larsen_alloc (double lambda1, double lambda2, const c_matrix *x, const c_vector *y);
void		larsen_free (larsen *l);

/* data.c */
double		larsen_centering_vector (c_vector *y);
c_vector	*larsen_centering_matrix (c_matrix *x);
c_vector	*larsen_normalizing_matrix (c_matrix *x);

/* larsen.c */
bool		larsen_regression_step (larsen *l);
bool		larsen_interpolate (larsen *l);
c_vector	*larsen_get_beta (larsen *l);
c_vector	*larsen_get_mu (larsen *l);
void		larsen_set_lambda1 (larsen *l, double t);
double		larsen_get_lambda1 (larsen *l, const bool scaled);

/* elasticnet.c */
bool		larsen_elasticnet (larsen *l, int maxiter);

/* bic.c */
double		larsen_eval_bic (larsen *l, double gamma);

#ifdef __cplusplus
}
#endif

#endif /* LARSEN_H_ */
