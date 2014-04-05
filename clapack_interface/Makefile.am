## Process this file with automake to produce Makefile.in

INCLUDES =\
	 -I./ -I../include

AM_CFLAGS =\
	 -Wall\
	 -g

lib_LTLIBRARIES = libclapack_interface.la

libclapack_interface_include_DATA = \
	cl_linalg.h\
	cl_matrix.h\
	cl_interface.h

EXTRA_DIST = $(liblarsen_include_DATA)

libclapack_interface_includedir = $(pkgincludedir)

libclapack_interface_la_SOURCES = \
	cl_matrix.c\
	cl_matrixops.c\
	cl_linalg.c

libclapack_interface_la_LDFLAGS = 

libclapack_interface_la_LIBADD = 

