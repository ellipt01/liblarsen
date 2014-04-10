## Process this file with automake to produce Makefile.in

INCLUDES =\
	 -I$(top_srcdir)/include -I$(top_srcdir)/clinalg -I$(top_srcdir)/.

AM_CFLAGS =\
	 -Wall\
	 -g

noinst_LTLIBRARIES = liblarsenexample.la 

noinst_HEADERS = example.h

liblarsenexample_la_SOURCES = \
	example.c\
	example_elasticnet.c

liblarsenexample_la_CPPFLAGS = -I$(top_srcdir)/.

liblarsenexample_la_LDFLAGS = -L$(top_srcdir)/clinalg -L$(top_srcdir)/src

liblarsenexample_la_LIBADD = -lclinalg -llarsen

