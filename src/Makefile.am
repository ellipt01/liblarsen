## Process this file with automake to produce Makefile.in

INCLUDES =\
	 -I$(top_srcdir)/include -I$(top_srcdir)/c_linalg

AM_CFLAGS =\
	 -Wall\
	 -g

lib_LTLIBRARIES = liblarsen.la

liblarsen_la_SOURCES = \
	activeset.c\
	data.c\
	equiangular.c\
	larsen.c\
	stepsize.c\
	utils.c\
	elasticnet.c\
	bic.c

liblarsen_la_CPPFLAGS = -I$(top_srcdir)/include

liblarsen_la_LDFLAGS = -L$(top_srcdir)/c_linalg

liblarsen_la_LIBADD = -lclinalg

