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

	size_t					n;	// number of data
	size_t					p;	// number of variables

	/* if true, loop of lasso or elastic net regression is terminated */
	bool					stop_loop;

	/* threshold for L1 penalty */
	double					lambda1;
	/* threshold for L2 penalty */
	double					lambda2;

	/* true: elastic net, false: lasso */
	bool					is_elnet;
	/* scale  = (is_elnet) ? 1 / sqrt(1 + lambda2) : 1
	 * scale2 = l->scale^2 */
	double					scale;
	double					scale2;

	const double			*y;		// data
	const double			*x;		// variables

	/* correlation */
	double					sup_c;	// sup(abs(c))
	double					*c;		// colleration: = X' * (y - mu)

	/* active set */
	activeset_oper		oper;
	size_t					sizeA;
	int						*A;		// active set
	int						*Ac;	// implement of A

	/* equiangular */
	double					absA;
	double					*u;		// equiangular vector
	double					*w;

	/* step_size */
	double					stepsize;	// step size which beta will be progress.

	/* solution */
	double					*beta;	// regression coefficients
	double					*mu;	// estimated response

	/* previous beta and mu */
	double					*beta_prev;	// previous beta
	double					*mu_prev;		// previous mu

	/* interpolation */
	bool					interp;		// interpolation was done or not.
	double					stepsize_intr;
	double					*beta_intr;	// interpolated beta
	double					*mu_intr;		// interpolated mu

	/* cholesky factorization */
	double					*chol;	// = chol(XA' * XA), where XA = X(A).

};

/* util.c */
larsen		*larsen_alloc (size_t n, size_t p, const double *y, const double *x, double lambda1, double lambda2);
void		larsen_free (larsen *l);

/* data.c */
double		larsen_centering_vector (c_vector *y);
c_vector	*larsen_centering_matrix (c_matrix *x);
c_vector	*larsen_normalizing_matrix (c_matrix *x);

/* larsen.c */
bool		larsen_regression_step (larsen *l);
bool		larsen_interpolate (larsen *l);
double		*larsen_get_beta (larsen *l);
double		*larsen_get_mu (larsen *l);
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
