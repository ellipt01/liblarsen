/*
 * larsen_private.h
 *
 *  Created on: 2014/05/14
 *      Author: utsugi
 */

#ifndef LARSEN_PRIVATE_H_
#define LARSEN_PRIVATE_H_

/* Private macros and entry of constants
 * which are only used internally */

/* cast pointer type to const int * */
#ifndef CINTP
#define CINTP(a)	((const int *) &(a))
#endif

/* min(a, b) */
#define LARSEN_MIN(a, b)	(((a) > (b)) ? (b) : (a))

/* posinf = +1. / 0. */
#ifndef LARSEN_POSINF
#define LARSEN_POSINF	(+1. / 0.)	// + infinity
#endif

/* following constants are set in larsen_linalg.c */
extern const int		ione;	//  1
extern const double	dzero;	//  0.
extern const double	dmone;	// -1.

#endif /* LARSEN_PRIVATE_H_ */
