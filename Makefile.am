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

# test and sample program
bin_PROGRAMS =\
	larsen-example\
	test_cholesky_decomp\
	test_cholesky_svx\
	test_cholesky_insert\
	test_cholesky_delete

# sample program
larsen_example_SOURCES = example/main.c
larsen_example_CPPFLAGS =\
	-I$(top_srcdir)/linalg\
	-I$(top_srcdir)/include\
	-I$(top_srcdir)/example
larsen_example_LDADD =\
	$(top_srcdir)/example/liblarsenexample.la

# test_cholesky_decomp
test_cholesky_decomp_SOURCES = test/test_cholesky.c test/test_cholesky_decomp.c
test_cholesky_decomp_CPPFLAGS =\
	-I$(top_srcdir)/linalg\
	-I$(top_srcdir)/include\
	-I$(top_srcdir)/test
test_cholesky_decomp_LDADD =\
	$(top_srcdir)/src/liblarsen.la

# test_cholesky_svx
test_cholesky_svx_SOURCES = test/test_cholesky.c test/test_cholesky_svx.c
test_cholesky_svx_CPPFLAGS =\
	-I$(top_srcdir)/linalg\
	-I$(top_srcdir)/include\
	-I$(top_srcdir)/test
test_cholesky_svx_LDADD =\
	$(top_srcdir)/src/liblarsen.la

# test_cholesky_insert
test_cholesky_insert_SOURCES = test/test_cholesky.c test/test_cholesky_insert.c
test_cholesky_insert_CPPFLAGS =\
	-I$(top_srcdir)/linalg\
	-I$(top_srcdir)/include\
	-I$(top_srcdir)/test
test_cholesky_insert_LDADD =\
	$(top_srcdir)/src/liblarsen.la

# test_cholesky_delete
test_cholesky_delete_SOURCES = test/test_cholesky.c test/test_cholesky_delete.c
test_cholesky_delete_CPPFLAGS =\
	-I$(top_srcdir)/linalg\
	-I$(top_srcdir)/include\
	-I$(top_srcdir)/test
test_cholesky_delete_LDADD =\
	$(top_srcdir)/src/liblarsen.la



# Copy all the spec files. Of cource, only one is actually used.
dist-hook:
	for specfile in *.spec; do \
		if test -f $$specfile; then \
			cp -p $$specfile $(distdir); \
		fi \
	done

