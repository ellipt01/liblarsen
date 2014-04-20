## Process this file with automake to produce Makefile.in

SUBDIRS = linalg include src example

liblarsendocdir = ${prefix}/doc/larsen
liblarsendoc_DATA = \
	README\
	COPYING\
	AUTHORS\
	ChangeLog\
	INSTALL\
	NEWS\
	TODO

EXTRA_DIST = $(liblarsendoc_DATA)

# Sample program
bin_PROGRAMS = larsen-example
larsen_example_SOURCES = example/main.c
larsen_example_CPPFLAGS =\
	-I$(top_srcdir)/linalg\
	-I$(top_srcdir)/include\
	-I$(top_srcdir)/example
larsen_example_LDADD =\
	$(top_srcdir)/example/liblarsenexample.la

# Copy all the spec files. Of cource, only one is actually used.
dist-hook:
	for specfile in *.spec; do \
		if test -f $$specfile; then \
			cp -p $$specfile $(distdir); \
		fi \
	done

