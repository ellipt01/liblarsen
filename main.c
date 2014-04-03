/*
 * main.c
 *
 *  Created on: 2014/03/17
 *      Author: utsugi
 */

#include <string.h>

#include "larsen.h"
#include "examples/examples.h"

char	fn[80] = "\0";
size_t	skipheaders = 0;
double	lambda2 = 0.;
double	start = 0.;
double	stop = 10.;
double	dt = 0.1;
int		maxsteps = 100;

void
usage (char *toolname)
{
	char	*p = strrchr (toolname, '/');
	if (p) p++;
	else p = toolname;

	fprintf (stderr, "\nUSAGE:\n%s -f <input_file>{:skipheaders}\n", p);
	fprintf (stderr, "\toptional  {-l <lambda2> -t <start>:<dt>:<stop> -m <maxsteps>}\n\n");
	exit (1);
}

bool
read_params (int argc, char **argv)
{
	int		i;
	bool	status = true;
	double	_start, _dt, _stop;

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

				case 'm':
					maxsteps = atoi (argv[++i]);
				break;

				default:
				break;
			}
		}
	}
	if (strlen (fn) <= 1) {
		fprintf (stderr, "ERROR: input file name is not specified\n");
		status = false;
	}
	if (_start >= _stop || floor ((_stop - _start) / _dt) <= 0) {
		fprintf (stderr, "ERROR: range of lambda1 invalid : %.2f:%.2f:%.2f\n", _start, _dt, _stop);
		status = false;
	}

	start = _start;
	dt = _dt;
	stop = _stop;

	return status;
}

void
fprintf_params (void)
{
	fprintf (stderr, "###########################################################\n\n");
	fprintf (stderr, "read file: \t\"%s\" (skip headers = %d)\n", fn, (int) skipheaders);
	fprintf (stderr, "lambda1 :\t[%.2f : %.2f : %.2f]\n", start, dt, stop);
	fprintf (stderr, "lambda2 :\t%.2f\n", lambda2);
	fprintf (stderr, "maxsteps :\t%d\n", maxsteps);
	fprintf (stderr, "\n###########################################################\n");
	return;
}

int
main (int argc, char **argv)
{
	larsen_data	*data = NULL;

	if (!read_params (argc, argv)) usage (argv[0]);
	fprintf_params ();

	data = read_data (fn, skipheaders);

	larsen_centering_vector (data->y);
	larsen_centering_matrix (data->x);
	larsen_normalizing_matrix (data->x);

	example_elasticnet (data->x, data->y, start, dt, stop, lambda2, maxsteps);

	larsen_data_free (data);

	return EXIT_SUCCESS;
}
