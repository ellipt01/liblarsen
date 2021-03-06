/*
 * main.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <larsen.h>
#include "example.h"

char		fn[80] = "\0";
size_t		skipheaders = 0;
double		lambda2 = 0.;
double		start = 0.;
double		stop = 100.;
double		dt = 0.1;
double		gamma_bic = 0.;	// traditional BIC
int			maxiter = 100;

void
usage (char *toolname)
{
	char	*p = strrchr (toolname, '/');
	if (p) p++;
	else p = toolname;

	fprintf (stderr, "\nUSAGE:\n%s -f <input_file>{:skipheaders} -l <lambda2> \n", p);
	fprintf (stderr, "[optional] { -t <start>:<dt>:<stop> -g <gamma of EBIC in [0, 1]> -m <maxsteps> }\n\n");
	exit (1);
}

bool
read_params (int argc, char **argv)
{
	int		i;
	bool	status = true;
	double	_start = start;
	double	_dt = dt;
	double	_stop = stop;
	double	_gamma = gamma_bic;

	if (argc <= 1) return false;

	for (i = 1; i < argc; i++) {
		char	*p;

		if (argv[i][0] == '-') {

			switch (argv[i][1]) {

				case 'f':
					p = strrchr (argv[++i], ':');
					if (p) {
						strcpy (fn, argv[i]);
						fn[strlen (argv[i]) - strlen (p)] = '\0';
						skipheaders = atoi (++p);
					} else strcpy (fn, argv[i]);
				break;

				case 'l':
					lambda2 = (double) atof (argv[++i]);
				break;

				case 't':
					sscanf (argv[++i], "%lf:%lf:%lf", &_start, &_dt, &_stop);
				break;

				case 'g':
					_gamma = (double) atof (argv[++i]);
				break;

				case 'm':
					maxiter = atoi (argv[++i]);
				break;

				default:
				break;
			}
		}
	}
	if (strlen (fn) <= 1) {
		fprintf (stderr, "ERROR: input file name is not specified.\n");
		status = false;
	}
	if (_start >= _stop || floor ((_stop - _start) / _dt) <= 0) {
		fprintf (stderr, "ERROR: range of lambda1 invalid : %.2f:%.2f:%.2f\n", _start, _dt, _stop);
		status = false;
	}
	if (_gamma < 0. || 1. < _gamma) {
		fprintf (stderr, "ERROR: gamma (%f) must be [0, 1].\n", _gamma);
		status = false;
	}

	start = _start;
	dt = _dt;
	stop = _stop;
	gamma_bic = _gamma;

	return status;
}

void
fprintf_params (void)
{
	fprintf (stderr, "###########################################################\n\n");
	fprintf (stderr, "read file: \t\"%s\" (skip headers = %d)\n", fn, (int) skipheaders);
	fprintf (stderr, "lambda1 :\t[%.2f : %.2f : %.2f]\n", start, dt, stop);
	fprintf (stderr, "lambda2 :\t%.2f\n", lambda2);
	fprintf (stderr, "maxiter :\t%d\n", maxiter);
	fprintf (stderr, "\n###########################################################\n");
	return;
}

int
main (int argc, char **argv)
{
	size_t		n;
	size_t		p;
	double		*y;
	double		*x;

	if (!read_params (argc, argv)) usage (argv[0]);
	fprintf_params ();

	read_data (fn, skipheaders, &n, &p, &y, &x);

	larsen_centering (n, 1, y);
	larsen_centering (n, p, x);
	larsen_normalizing (n, p, x);

	example_elasticnet (n, p, x, y, start, dt, stop, lambda2, gamma_bic, maxiter);

	free (y);
	free (x);

	return EXIT_SUCCESS;
}
