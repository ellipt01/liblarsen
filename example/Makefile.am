## Process this file with automake to produce Makefile.in

# sample program
noinst_PROGRAMS = larsen-example

noinst_LTLIBRARIES = liblarsenexample.la 

noinst_HEADERS = example.h

liblarsenexample_la_SOURCES =\
	example.c\
	example_elasticnet.c

liblarsenexample_la_CPPFLAGS =\
	-I$(top_srcdir)/include\
	-I$(top_srcdir)/example

liblarsenexample_la_LDFLAGS =\
	-L$(top_srcdir)/src\
	-L$(top_srcdir)/src/linreg

liblarsenexample_la_LIBADD = -llinreg -llarsen

# sample program
larsen_example_SOURCES = main.c
larsen_example_CPPFLAGS =\
	-I$(top_srcdir)/include\
	-I$(top_srcdir)/example

larsen_example_LDADD =\
	$(top_srcdir)/example/liblarsenexample.la
