## Process this file with automake to produce Makefile.in

SUBDIRS = src/linsys src include example test

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

# Copy all the spec files. Of cource, only one is actually used.
dist-hook:
	for specfile in *.spec; do \
		if test -f $$specfile; then \
			cp -p $$specfile $(distdir); \
		fi \
	done

