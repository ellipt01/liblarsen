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
#include <cl_interface.h>

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

	// if true, loop of larsen regression will be terminated.
	bool					stop_loop;

	/* threshold for L1 penalty */
	double					lambda1;
	/* threshold for L2 penalty */
	double					lambda2;

	// false: lasso, and true: elastic net estimator.
	bool					do_scaling;
	// = 1 / sqrt(1 + lambda2)
	double					scale;

	const cl_matrix		*x;		// predictors
	const cl_vector		*y;		// response

	/* correlation */
	double					sup_c;	// sup(abs(c))
	cl_vector				*c;		// colleration: = X' * (y - mu)

	/* active set */
	activeset_oper		oper;
	cl_vector_int			*A;		// active set
	cl_vector_int			*Ac;	// implement of A

	/* equiangular */
	double					absA;
	cl_vector				*u;		// equiangular vector
	cl_vector				*w;

	/* step_size */
	double					stepsize;	// step size which beta will be progress.

	/* solution */
	cl_vector				*beta;	// regression coefficients
	cl_vector				*mu;	// estimated response

	/* backup of solution */
	cl_vector				*beta_prev;	// previous beta
	cl_vector				*mu_prev;		// previous mu

	/* interpolation */
	bool					interp;		// interpolation was done or not.
	double					stepsize_intr;
	cl_vector				*beta_intr;	// interpolated beta
	cl_vector				*mu_intr;		// interpolated mu

	/* cholesky factorization */
	cl_matrix				*chol;	// = chol(XA' * XA), where XA = X(A).

};

/* util.c */
larsen		*larsen_alloc (double lambda1, double lambda2, const cl_matrix *x, const cl_vector *y);
void		larsen_free (larsen *l);

/* data.c */
double		larsen_centering_vector (cl_vector *y);
cl_vector	*larsen_centering_matrix (cl_matrix *x);
cl_vector	*larsen_normalizing_matrix (cl_matrix *x);

/* larsen.c */
bool		larsen_regression_step (larsen *l);
bool		larsen_interpolate (larsen *l);
cl_vector	*larsen_get_beta (larsen *l);
cl_vector	*larsen_get_mu (larsen *l);
double		larsen_increment_lambda1 (larsen *l, double dt);

/* elasticnet.c */
bool		larsen_elasticnet (larsen *l, int maxiter);

/* bic.c */
double		larsen_eval_bic (larsen *l, double gamma);

#ifdef __cplusplus
}
#endif

#endif /* LARSEN_H_ */
