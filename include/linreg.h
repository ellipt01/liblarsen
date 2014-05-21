/*
 * linreg.h
 *
 *  Wrapper of lapack and qrupdate
 *
 *  Created on: 2014/04/08
 *      Author: utsugi
 */

#ifndef LINREG_H_
#define LINREG_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdbool.h>

#ifndef LINREG_INDEX_OF_MATRIX
#define LINREG_INDEX_OF_MATRIX(i, j, lda) ((i) + (j) * (lda))
#endif

typedef struct s_linreg		linreg;
typedef struct s_penalty		penalty;

/* structure of linear regression model
 *   b = Z * beta,
 *   b = [y; 0]
 *   Z = scale * [X; sqrt(lambda2) * J]
 */
struct s_linreg {
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

	/* penalty term.
	 * if pen == NULL && lambda2 > 0, ridge regression is assumed. */
	const penalty			*pen;

};

/* penalty term */
struct s_penalty {
	size_t					pj;		// rows of r
	size_t					p;		// columns of r
	double					*r;		// pj x p penalty matrix
};

/* linreg.c */
linreg			*linreg_alloc (const double lambda2, const size_t n, const size_t p, const double *y, const double *x);
void			linreg_free (linreg *l);

void			linreg_centering_y (linreg *lreg);
void			linreg_centering_x (linreg *lreg);
void			linreg_normalizing_x (linreg *lreg);
void			linreg_standardizing_x (linreg *lreg);

penalty		*penalty_alloc (const size_t p1, const size_t p, const double *r);
void			penalty_free (penalty *pen);

void			linreg_set_penalty (linreg *lreg, const double a, const double b, const penalty *pen);

#ifdef __cplusplus
}
#endif

#endif /* LINREG_H_ */
