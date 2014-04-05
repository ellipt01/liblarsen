/*
 * example.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <string.h>

#include <larsen.h>
#include "example.h"

larsen_data *
larsen_data_alloc (const size_t size1, const size_t size2)
{
	larsen_data	*data = (larsen_data *) malloc (sizeof (larsen_data));
 	data->x = cl_matrix_alloc (size1, size2);
	data->y = cl_vector_alloc (size1);
	return data;
}

void
larsen_data_free (larsen_data *data)
{
	if (data) {
		if (!cl_matrix_is_empty (data->x)) cl_matrix_free (data->x);
		if (!cl_vector_is_empty (data->y)) cl_vector_free (data->y);
		free (data);
	}
	return;
}

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

larsen_data *
read_data (char *fn, int skip_header)
{
	int			i, j;
	size_t		size1;
	size_t		size2;
	larsen_data	*data = NULL;

	char		buf[100 * BUFSIZ];
	FILE		*fp;

	count_data (fn, skip_header, &size1, &size2);

	data = larsen_data_alloc (size1, size2);

	if ((fp = fopen (fn, "r")) == NULL) return NULL;
	i = 0;
	while (fgets (buf, 100 * BUFSIZ, fp) != NULL) {
		char	*p;
		if (buf[0] == '#' || buf[0] == '\n') continue;
		if (i - skip_header >= 0) {
			for (j = 0, p = strtok (buf, "\t "); p != NULL; j++, p = strtok (NULL, "\t ")) {
				double	val = (double) atof (p);
				if (j >= data->x->size2) cl_vector_set (data->y, i - skip_header, val);
				else cl_matrix_set (data->x, i - skip_header, j, val);
			}
		}
		i++;
	}
	fclose (fp);
	return data;
}

void
output_solutionpath (int iter, larsen *l)
{
	int			i;
	char		fn[80];
	FILE		*fp;
	cl_vector	*beta = larsen_get_beta (l);

	for (i = 0; i < beta->size; i++) {

		sprintf (fn, "beta%03d.res", i);

		if (iter == 0) fp = fopen (fn, "w");
		else fp = fopen (fn, "aw");
		if (fp == NULL) continue;

		fprintf (fp, "%d\t%.4e\t%.4e\n", iter, l->lambda1, cl_vector_get (beta, i));
		fclose (fp);
	}
	cl_vector_free (beta);
	return;
}

void
fprintf_beta (FILE *stream, int iter, larsen *l)
{
	cl_vector	*beta = larsen_get_beta (l);

	fprintf (stream, "beta[%04d, t = %.2f] = \n", iter, l->lambda1);
	cl_vector_fprintf (stream, beta, "%.4e");
	fprintf (stream, "]\n");

	cl_vector_free (beta);

	return;
}

