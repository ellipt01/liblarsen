/*
 * larsen_private.h
 *
 *  Created on: 2014/05/14
 *      Author: utsugi
 */

#ifndef LARSEN_PRIVATE_H_
#define LARSEN_PRIVATE_H_

#ifndef CINTP
#define CINTP(a)	((const int *) &(a))
#endif

#define LARSEN_MIN(a, b)	(((a) > (b)) ? (b) : (a))

#ifndef LARSEN_POSINF
#define LARSEN_POSINF	(1. / 0.)	// + infinity
#endif

extern const int		ione;
extern const double	dzero;
extern const double	dmone;

#endif /* LARSEN_PRIVATE_H_ */
