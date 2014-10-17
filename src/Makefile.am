## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

AM_CFLAGS =\
	 -Wall\
	 -g

lib_LTLIBRARIES = liblarsen.la

liblarsen_la_SOURCES = \
	private.c\
	atomic.c\
	mmio.c\
	mmreal.c\
	linregmodel.c\
	linalg.c\
	activeset.c\
	equiangular.c\
	larsen.c\
	stepsize.c\
	update.c\
	elasticnet.c\
	bic.c

liblarsen_la_CPPFLAGS = -I$(top_srcdir)/include

liblarsen_la_LDFLAGS =

liblarsen_la_LIBADD =

