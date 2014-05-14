/*
 * larsen_linalg.c
 *
 *  Created on: 2014/05/04
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <larsen_linalg.h>

#include "larsen_private.h"

const int		ione  =  1;
const double	dzero =  0.;
const double	dmone = -1.;

void
larsen_linalg_error (const char * function_name, const char *error_msg)
{
	fprintf (stderr, "ERROR: %s: %s\n", function_name, error_msg);
	exit (1);
}

void
matrix_add_rowcols (const size_t m, const size_t n, double **a, const size_t dm, const size_t dn)
{
	int		j;
	if (dm <= 0 && dn <= 0) return;
	if (*a == NULL) {
		*a = (double *) malloc ((m + dm) * (n + dn) * sizeof (double));
		return;
	}

	*a = (double *) realloc (*a, (m + dm) * (n + dn) * sizeof (double));

	if (dm > 0) {
		double	*col = (double *) malloc (m * sizeof (double));
		for (j = n; 0 < j; j--) {
			dcopy_ (CINTP (m), *a + LARSEN_INDEX_OF_MATRIX (0, j, m), &ione, col, &ione);
			dcopy_ (CINTP (m), col, &ione, *a + LARSEN_INDEX_OF_MATRIX (0, j, m + dm), &ione);
		}
		free (col);
	}
	return;
}

void
matrix_remove_rowcols (const size_t m, const size_t n, double **a, const size_t dm, const size_t dn)
{
	int		j;
	if (dm <= 0 && dn <= 0) return;
	if (m - dm <= 0 || n - dn <= 0) {
		if (*a) free (*a);
		*a = NULL;
		return;
	}
	if (dm > 0) {
		double	*col = (double *) malloc ((m - dm) * sizeof (double));
		for (j = 1; j < n - dn; j++) {
			int		mm = (int) (m - dm);
			dcopy_ (CINTP (mm), *a + LARSEN_INDEX_OF_MATRIX (0, j, m), &ione, col, &ione);
			dcopy_ (CINTP (mm), col, &ione, *a + LARSEN_INDEX_OF_MATRIX (0, j, m - dm), &ione);
		}
		free (col);
	}

	*a = (double *) realloc (*a, (m - dm) * (n - dn) * sizeof (double));

	return;
}
