## Process this file with automake to produce Makefile.in

INCLUDES =\
	 -I./ -I../include

AM_CFLAGS =\
	 -Wall\
	 -g

lib_LTLIBRARIES = libclinalg.la

libclinalg_include_DATA = \
	c_matrix.h\
	c_linalg_cholesky.h\
	c_linalg.h

EXTRA_DIST = $(liblarsen_include_DATA)

libclinalg_includedir = $(pkgincludedir)

libclinalg_la_SOURCES = \
	c_matrix.c\
	c_matrixops.c\
	c_linalg_cholesky.c

libclinalg_la_LDFLAGS = 

libclinalg_la_LIBADD = 

