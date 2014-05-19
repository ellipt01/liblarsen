/*
 * linsys_private.h
 *
 *  Created on: 2014/05/14
 *      Author: utsugi
 */

#ifndef LINSYS_PRIVATE_H_
#define LINSYS_PRIVATE_H_

/* Private macros and entry of constants
 * which are only used internally */

/* cast pointer type to const int * */
#ifndef LINSYS_CINTP
#define LINSYS_CINTP(a)	((const int *) &(a))
#endif

/* min(a, b) */
#define LINSYS_MIN(a, b)	(((a) > (b)) ? (b) : (a))

/* posinf = +1. / 0. */
#ifndef LINSYS_POSINF
#define LINSYS_POSINF	(+1. / 0.)	// + infinity
#endif

/* following constants are set in linsys.c */
extern const int		ione;	//  1
extern const double	dzero;	//  0.
extern const double	dmone;	// -1.

void	linsys_error (const char * function_name, const char *error_msg, const char *file, const int line);

#endif /* LINSYS_PRIVATE_H_ */
