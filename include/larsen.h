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

#include <linreg.h>

typedef enum {
	ACTIVESET_ACTION_NONE	= -1,
	ACTIVESET_ACTION_ADD		=  0,	// add new variable to the active set
	ACTIVESET_ACTION_DROP	=  1,	// remove a variable from the active set
	NUM_ACTIVESET_ACTION		=  2
} ActiveSetAction;

typedef struct {
	ActiveSetAction		action;
	int						index_of_A;	// position of A
	int						column_of_X;	// operand column of matrix X
} activeset_operation;

typedef struct s_larsen	larsen;

struct s_larsen {

	/* if true, loop of regression is terminated */
	bool					stop_loop;

	/* linear system of regression equations */
	const linreg			*lreg;

	/* threshold for L1 penalty */
	double					lambda1;

	/* correlation */
	double					sup_c;	// sign : sup(abs(c))
	double					*c;		// correlation : = X' * (y - mu)

	/* active set */
	activeset_operation	oper;
	size_t					sizeA;
	int						*A;		// active set

	/* equiangular */
	double					absA;
	double					*u;		// equiangular vector
	double					*w;

	/* step_size */
	double					stepsize;	// step size which beta will be progress

	/* solution */
	double					nrm1;	// sum_j |beta(j)|
	double					*beta;	// regression coefficients
	double					*mu;	// estimation of the response y

	/* previous beta and mu */
	double					nrm1_prev;		// sum_j |beta_prev(j)|
	double					*beta_prev;	// backup of previous beta
	double					*mu_prev;		// backup of previous mu

	/* cholesky factorization */
	double					*chol;	// = chol(XA' * XA), where XA = X(A)

};

/* linalg.c */
int			larsen_linalg_cholesky_svx (const size_t size, double *l, const size_t lda, double *b);
int			larsen_linalg_cholesky_insert (const size_t n, double **r, const int index, double *u);
void		larsen_linalg_cholesky_delete (const size_t size, double **r, const int index);

/* util.c */
larsen		*larsen_alloc (const linreg *lreg, const double lambda1);
void		larsen_free (larsen *l);
double		*larsen_copy_beta (const larsen *l, bool scaling);
double		*larsen_copy_mu (const larsen *l, bool scaling);
void		larsen_set_lambda1 (larsen *l, double t);
double		larsen_get_lambda1 (const larsen *l, bool scaling);
bool		larsen_is_regtype_lasso (const larsen *l);
bool		larsen_is_regtype_ridge (const larsen *l);

/* data.c */
double		*larsen_centering (const size_t size1, const size_t size2, double *x);
double		*larsen_normalizing (const size_t size1, const size_t size2, double *x);

/* larsen.c */
bool		larsen_regression_step (larsen *l);
bool		larsen_does_need_interpolation (const larsen *l);

/* elasticnet.c */
bool		larsen_elasticnet (larsen *l, int maxiter);

/* bic.c */
double		larsen_eval_bic (const larsen *l, double gamma);

#ifdef __cplusplus
}
#endif

#endif /* LARSEN_H_ */
