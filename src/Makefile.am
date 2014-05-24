## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

AM_CFLAGS =\
	 -Wall\
	 -g

lib_LTLIBRARIES = liblarsen.la

liblarsen_la_SOURCES = \
	linalg.c\
	activeset.c\
	equiangular.c\
	larsen.c\
	stepsize.c\
	utils.c\
	estimator.c\
	bic.c

liblarsen_la_CPPFLAGS = -I$(top_srcdir)/include

liblarsen_la_LDFLAGS = -L$(top_srcdir)/src/linreg

liblarsen_la_LIBADD = -llinreg

