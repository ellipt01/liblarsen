## Process this file with automake to produce Makefile.in

noinst_PROGRAMS =\
	test_cholesky_svx\
	test_cholesky_insert\
	test_cholesky_delete

noinst_HEADERS = test_larsen.h

# test_cholesky_svx
test_cholesky_svx_SOURCES = test_cholesky.c test_cholesky_svx.c
test_cholesky_svx_CPPFLAGS =\
	-I$(top_srcdir)/include\
	-I$(top_srcdir)/test
test_cholesky_svx_LDADD =\
	$(top_srcdir)/src/larsen/liblarsen.la

# test_cholesky_insert
test_cholesky_insert_SOURCES = test_cholesky.c test_cholesky_insert.c
test_cholesky_insert_CPPFLAGS =\
	-I$(top_srcdir)/include\
	-I$(top_srcdir)/test
test_cholesky_insert_LDADD =\
	$(top_srcdir)/src/larsen/liblarsen.la

# test_cholesky_delete
test_cholesky_delete_SOURCES = test_cholesky.c test_cholesky_delete.c
test_cholesky_delete_CPPFLAGS =\
	-I$(top_srcdir)/include\
	-I$(top_srcdir)/test
test_cholesky_delete_LDADD =\
	$(top_srcdir)/src/larsen/liblarsen.la
