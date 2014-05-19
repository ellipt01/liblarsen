#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([liblarsen], [0.1])
AC_CONFIG_SRCDIR([example/example.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_LIBTOOL

AM_INIT_AUTOMAKE

# Checks for libraries.
AC_CHECK_LIB([m], [floor], [])
AC_CHECK_LIB([gfortran], [main],[])

AC_CHECK_LIB(
	[qrupdate],
	[dchinx_],
	[],
	[AC_MSG_ERROR("ERROR: cannot find library for -lqrupdate. check LDFLAGS and retry - stop.")]
)
AC_CHECK_FUNCS([dchinx_ dchdex_])


# Checks for header files.
AC_CHECK_HEADER([stdlib.h])
AC_CHECK_HEADER([string.h])
AC_CHECK_HEADER([lapack.h])
AC_CHECK_HEADER([blas.h])
AC_CHECK_HEADER([qrupdate.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor pow sqrt log strrchr])
AC_CHECK_FUNCS([dpotrf_ dpotrs_])

# pkgconfig
LIBLARSEN_DESCRIPTION='C implementation of LARS-EN algorithm'
AC_SUBST([LIBLARSEN_DESCRIPTION])
AC_CONFIG_FILES([liblarsen.pc:liblarsen.pc.in])

AC_CONFIG_FILES([Makefile
                 include/Makefile
                 src/linsys/Makefile
                 src/Makefile
                 example/Makefile
                 test/Makefile])
AC_OUTPUT
