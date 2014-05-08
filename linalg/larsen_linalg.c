/*
 * larsen_linalg.c
 *
 *  Created on: 2014/05/04
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <larsen_linalg.h>

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
			cblas_dcopy (m, *a + LARSEN_INDEX_OF_MATRIX (0, j, m), 1, col, 1);
			cblas_dcopy (m, col, 1, *a + LARSEN_INDEX_OF_MATRIX (0, j, m + dm), 1);
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
			cblas_dcopy (m - dm, *a + LARSEN_INDEX_OF_MATRIX (0, j, m), 1, col, 1);
			cblas_dcopy (m - dm, col, 1, *a + LARSEN_INDEX_OF_MATRIX (0, j, m - dm), 1);
		}
		free (col);
	}

	*a = (double *) realloc (*a, (m - dm) * (n - dn) * sizeof (double));

	return;
}
