## Process this file with automake to produce Makefile.in

SUBDIRS = c_linalg include src example

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

bin_PROGRAMS = larsen-example
larsen_example_SOURCES = example/main.c
larsen_example_CPPFLAGS = -I$(top_srcdir)/c_linalg -I$(top_srcdir)/include -I.
larsen_example_LDADD = \
	$(top_srcdir)/c_linalg/libclinalg.la\
	$(top_srcdir)/src/liblarsen.la\
	$(top_srcdir)/example/liblarsenexample.la -lqrupdate -lm

# Copy all the spec files. Of cource, only one is actually used.
dist-hook:
	for specfile in *.spec; do \
		if test -f $$specfile; then \
			cp -p $$specfile $(distdir); \
		fi \
	done

