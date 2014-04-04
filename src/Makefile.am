## Process this file with automake to produce Makefile.in

INCLUDES =\
	 -I../include -I../clapack_interface

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

liblarsen_la_LDFLAGS = 

liblarsen_la_LIBADD = \
	 -lclapack_interface

