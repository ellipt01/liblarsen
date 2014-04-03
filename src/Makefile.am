## Process this file with automake to produce Makefile.in

## Created by Anjuta - will be overwritten
## If you don't want it to overwrite it,
## 	Please disable it in the Anjuta project configuration

INCLUDES =\
	 -I../include -I/usr/local/include/clapack_interface -I../include

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

