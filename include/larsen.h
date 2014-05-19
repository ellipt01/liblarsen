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

#include <linsys.h>

/* lapack */
#ifdef HAVE_LAPACK_H
#include <lapack.h>
#else
extern double	dlamch_ (const char *cmach);
extern void	dpotrs_ (const char *uplo, const int *n, const int *nrhs, const double *a, const int *lda, double *b, const int *ldb, int *info);
#endif

/* qrupdate: choesky linsert/delete */
#ifdef HAVE_QRUPDATE_H
#include <qrupdate.h>
#else
extern void	dchinx_ (const int *n, double *R, const int *ldr, const int *j, const double *u, double *w, int *info);
extern void	dchdex_ (const int *n, double *R, const int *ldr, const int *j, double *w);
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

	/* linear system of regression equations */
	const linsys			*sys;

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
	double					*beta;	// regression coefficients
	double					*mu;	// estimated response

	/* previous beta and mu */
	double					*beta_prev;	// previous beta
	double					*mu_prev;		// previous mu

	/* interpolation */
	bool					interp;		// interpolation was done or not
	double					stepsize_intr;
	double					*beta_intr;	// interpolated beta
	double					*mu_intr;		// interpolated mu

	/* cholesky factorization */
	double					*chol;	// = chol(XA' * XA), where XA = X(A)

};

/* linalg.c */
int			larsen_linalg_cholesky_svx (const size_t size, double *l, const size_t lda, double *b);
int			larsen_linalg_cholesky_insert (const size_t n, double **r, const int index, double *u);
void		larsen_linalg_cholesky_delete (const size_t size, double **r, const int index);

/* util.c */
larsen		*larsen_alloc (const linsys *sys, const double lambda1, const double lambda2);
void		larsen_free (larsen *l);
double		*larsen_copy_beta_navie (const larsen *l);
double		*larsen_copy_beta_elasticnet (const larsen *l);
double		*larsen_copy_mu_navie (const larsen *l);
double		*larsen_copy_mu_elasticnet (const larsen *l);
void		larsen_set_lambda1 (larsen *l, double t);

/* data.c */
double		*larsen_centering (const size_t size1, const size_t size2, double *x);
double		*larsen_normalizing (const size_t size1, const size_t size2, double *x);

/* larsen.c */
bool		larsen_regression_step (larsen *l);
bool		larsen_interpolate (larsen *l);

/* elasticnet.c */
bool		larsen_elasticnet (larsen *l, int maxiter);

/* bic.c */
double		larsen_eval_bic (const larsen *l, double gamma);

#ifdef __cplusplus
}
#endif

#endif /* LARSEN_H_ */
