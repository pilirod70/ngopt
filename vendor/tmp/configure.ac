#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(idba, 0.19, ypeng@cs.hku.hk)
AC_CONFIG_SRCDIR([src/common/BitOperation.cpp])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE(idba, 0.19)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
#AC_CHECK_LIB(boost_filesystem, main)
#AC_CHECK_LIB(boost_program_options, main)
AC_PROG_RANLIB

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST

# Checks for library functions.
AC_FUNC_MEMCMP
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([sqrt])

AC_OUTPUT(Makefile bin/Makefile lib/Makefile test/Makefile)
