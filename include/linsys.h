/*
 * linsys.h
 *
 *  Wrapper of lapack and qrupdate
 *
 *  Created on: 2014/04/08
 *      Author: utsugi
 */

#ifndef LINSYS_H_
#define LINSYS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdbool.h>

#ifndef LINSYS_INDEX_OF_MATRIX
#define LINSYS_INDEX_OF_MATRIX(i, j, lda) ((i) + (j) * (lda))
#endif

typedef struct s_linsys		linsys;
typedef struct s_penalty		penalty;

typedef enum {
	REGULARIZATION_LASSO			= 0,
	REGULARIZATION_RIDGE			= 1,
	REGULARIZATION_USER_DEFINED	= 2
} RegularizationType;

/* structure for linear system of regression equations */
struct s_linsys {
	size_t					n;	// number of data
	size_t					p;	// number of variables

	double					*y;		// data
	double					*x;		// variables

	bool					ycentered;
	bool					xcentered;
	bool					xnormalized;

	double					*meany;	// mean(y)
	double					*meanx;	// meanx[j] = mean( X(:,j) )
	double					*normx;	// normx[j] = norm( X(:,j) )

	/* threshold for L2 penalty */
	double					lambda2;

	/* scale factor of system */
	double					scale2;	// scale2 = 1 / (a + b * lambda2)
	double					scale;		// scale = sqrt(scale2)

	/* use default regularization type (= ridge regression) */
	RegularizationType	regtype;
	/* penalty term.
	 * if pen == NULL && lambda2 > 0, ridge regression is assumed. */
	const penalty			*pen;

};

/* penalty term */
struct s_penalty {
	size_t					pj;		// rows of r
	size_t					p;		// cols of r
	double					*r;		// pj x p penalty matrix
};

/* linsys.c */
linsys			*linsys_alloc (const double lambda2, const size_t n, const size_t p, const double *y, const double *x);
void			linsys_free (linsys *l);

void			linsys_centering_y (linsys *lsys);
void			linsys_centering_x (linsys *lsys);
void			linsys_normalizing_x (linsys *lsys);
void			linsys_standardizing_x (linsys *lsys);

penalty		*penalty_alloc (const size_t p1, const size_t p, const double *r);
void			penalty_free (penalty *pen);

void			linsys_set_penalty (linsys *lsys, const double a, const double b, const penalty *pen);

bool			linsys_is_regtype_lasso (const linsys *lsys);
bool			linsys_is_regtype_ridge (const linsys *lsys);

double			linsys_get_lambda2 (const linsys *lsys);
double			linsys_get_scale (const linsys *lsys);
double			linsys_get_scale2 (const linsys *lsys);

const size_t	linsys_get_n (const linsys *lsys);
const size_t	linsys_get_p (const linsys *lsys);
const size_t	linsys_get_pj (const linsys *lsys);

const double	*linsys_get_x (const linsys *lsys);
const double	*linsys_get_y (const linsys *lsys);
const double	*linsys_get_penalty (const linsys *lsys);

#ifdef __cplusplus
}
#endif

#endif /* LINSYS_H_ */
