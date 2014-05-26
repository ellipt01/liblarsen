/*
 * example.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <larsen.h>

static void
count_data (char *fn, int skip_header, size_t *row, size_t *col)
{
	size_t	ndata = 0;
	size_t	npred = 0;
	char	buf[100 * BUFSIZ];
	FILE	*fp;

	*row = 0;
	*col = 0;

	if ((fp = fopen (fn, "r")) == NULL) {
		fprintf (stderr, "ERROR: cannot open file %s.\n", fn);
		exit (1);
	}
	while (fgets (buf, 100 * BUFSIZ, fp) != NULL) {
		if (buf[0] == '#' || buf[0] == '\n') continue;
		if (ndata - skip_header == 0) {
			char	*p;
			for (p = strtok (buf, "\t "); p != NULL; p = strtok (NULL, "\t ")) npred++;
		}
		ndata++;
	}
	fclose (fp);
	*row = ndata - skip_header;
	*col = npred - 1;
	return;
}

void
read_data (char *fn, int skip_header, size_t *n, size_t *p, double **y, double **x)
{
	int			i, j;
	size_t		size1;
	size_t		size2;

	char		buf[100 * BUFSIZ];
	FILE		*fp;

	double		*_y;
	double		*_x;

	count_data (fn, skip_header, &size1, &size2);

	_y = (double *) malloc (size1 * sizeof (double));
	_x = (double *) malloc (size1 * size2 * sizeof (double));

	if ((fp = fopen (fn, "r")) == NULL) return;
	i = 0;
	while (fgets (buf, 100 * BUFSIZ, fp) != NULL) {
		char	*p;
		if (buf[0] == '#' || buf[0] == '\n') continue;
		if (i - skip_header >= 0) {
			for (j = 0, p = strtok (buf, "\t "); p != NULL; j++, p = strtok (NULL, "\t ")) {
				double	val = (double) atof (p);
				if (j >= size2) _y[i - skip_header] = val;
				else _x[i - skip_header + j * size1] = val;
			}
		}
		i++;
	}
	fclose (fp);

	*n = size1;
	*p = size2;
	*x = _x;
	*y = _y;

	return;
}

void
fprintf_beta (FILE *stream, int iter, larsen *l)
{
	int		i;
	size_t	p = l->lreg->p;
	double	*beta = larsen_copy_beta (l, true);

	fprintf (stream, "beta[%04d, t = %.2f] = \n", iter, l->lambda1);
	for (i = 0; i < p; i++) fprintf (stream, "%.4e\n", beta[i]);
	fprintf (stream, "]\n");
	free (beta);

	return;
}

