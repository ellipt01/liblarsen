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

#include <stddef.h>
#include <stdbool.h>

#ifndef LARSEN_INDEX_OF_MATRIX
#define LARSEN_INDEX_OF_MATRIX(i, j, lda) ((i) + (j) * (lda))
#endif

typedef enum {
	ACTIVESET_ACTION_NONE	= -1,
	ACTIVESET_ACTION_ADD		=  0,	// add new variable to the active set
	ACTIVESET_ACTION_DROP	=  1,	// remove a variable from the active set
	NUM_ACTIVESET_ACTION		=  2
} ActiveSetAction;

typedef struct {
	ActiveSetAction	action;
	int					index_of_A;	// position of A
	int					column_of_X;	// operand column of matrix X
} activeset_operation;

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
	bool					is_lasso;

	/* scale  = (is_lasso) ? 1 : 1 / sqrt(1 + lambda2)
	 * scale2 = l->scale^2 */
	double					scale;
	double					scale2;

	const double			*y;		// data
	const double			*x;		// variables

	/* correlation */
	double					sup_c;	// sup(abs(c))
	double					*c;		// colleration: = X' * (y - mu)

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
	double					nrm1;
	double					*beta;	// regression coefficients
	double					*mu;	// estimated response

	/* previous beta and mu */
	double					nrm1_prev;
	double					*beta_prev;	// previous beta
	double					*mu_prev;		// previous mu

	/* cholesky factorization */
	double					*chol;	// = chol(XA' * XA), where XA = X(A)

};

/* util.c */
larsen		*larsen_alloc (const size_t n, const size_t p, const double *y, const double *x, const double lambda1, const double lambda2);
void		larsen_free (larsen *l);
double		*larsen_copy_beta (const larsen *l, bool scale);
double		*larsen_copy_mu (const larsen *l, bool scale);
void		larsen_set_lambda1 (larsen *l, double t);
double		larsen_get_lambda1 (const larsen *l, bool scaling);

/* data.c */
double		*larsen_centering (const size_t size1, const size_t size2, double *x);
double		*larsen_normalizing (const size_t size1, const size_t size2, double *x);

/* larsen.c */
bool		larsen_regression_step (larsen *l);
bool		larsen_does_need_interpolation (const larsen *l);
//bool		larsen_interpolate (larsen *l);

/* elasticnet.c */
bool		larsen_elasticnet (larsen *l, int maxiter);

/* bic.c */
double		larsen_eval_bic (const larsen *l, double gamma);

#ifdef __cplusplus
}
#endif

#endif /* LARSEN_H_ */
