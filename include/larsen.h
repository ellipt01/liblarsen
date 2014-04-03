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

#include <cl_interface.h>

typedef enum {
	ACTIVESET_ACTION_NONE	= -1,
	ACTIVESET_ACTION_ADD		=  0,
	ACTIVESET_ACTION_DROP	=  1,
	NUM_ACTIVESET_ACTION		=  2
} ActiveSetAction;

typedef struct {
	ActiveSetAction	action;
	int					column;	// index of operand column of matrix X
	int					index;		// position of added / removed item in A
} activeset_oper;

typedef struct s_larsen	larsen;

struct s_larsen {

	bool					stop_loop;

	/* threshold for L1 penalty */
	double					lambda1;
	/* threshold for L2 penalty */
	double					lambda2;

	bool					do_scaling;
	double					scale;

	const cl_matrix		*x;
	const cl_vector		*y;

	/* correlation */
	double					sup_c;
	cl_vector				*c;

	/* active set */
	activeset_oper		oper;
	cl_vector_int			*A;
	cl_vector_int			*Ac;

	/* equiangular */
	double					absA;
	cl_vector				*u;
	cl_vector				*w;

	/* step_size */
	double					stepsize;

	/* solution */
	cl_vector				*beta;
	cl_vector				*mu;

	/* backup of solution */
	cl_vector				*beta_prev;
	cl_vector				*mu_prev;

	/* interpolation */
	bool					interp;
	double					stepsize_intr;
	cl_vector				*beta_intr;
	cl_vector				*mu_intr;

	/* cholesky factorization */
	cl_matrix				*chol;

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
void		larsen_increment_lambda1 (larsen *l, double dt);
bool		larsen_loop_continue (larsen *l, double stop);

/* elasticnet.c */
bool		larsen_elasticnet (larsen *l, int maxiter);

/* bic.c */
double		larsen_eval_bic (larsen *l, double gamma);

#ifdef __cplusplus
}
#endif

#endif /* LARSEN_H_ */
