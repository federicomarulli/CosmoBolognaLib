AC_INIT([libsharp],[1.0.0])
AM_INIT_AUTOMAKE([foreign subdir-objects -Wall -Werror])
AM_MAINTAINER_MODE([enable])


dnl
dnl Needed for linking on Windows.
dnl Protect with m4_ifdef because AM_PROG_AR is required in
dnl autoconf >= 1.12 when using -Wall, but the macro is
dnl absent in old versions of autoconf.
dnl
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

LT_INIT
AC_CONFIG_MACRO_DIR([m4])

dnl
dnl Enable silent build rules if this version of Automake supports them
dnl
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])


m4_version_prereq(2.70, [AC_PROG_CC], [AC_PROG_CC_C99])
AC_OPENMP
CPPFLAGS="$CPPFLAGS $OPENMP_CFLAGS"
CFLAGS="$CFLAGS $OPENMP_CFLAGS"

# adding the lib to the files to link
LIBS="-lm"

LT_INIT

tmpval=`echo $CFLAGS | grep -c '\-DMULTIARCH'`
AM_CONDITIONAL([HAVE_MULTIARCH], [test $tmpval -gt 0])

PACKAGE_LIBS="-lsharp"

dnl
dnl Create pkgconfig .pc file.
dnl
AX_CREATE_PKGCONFIG_INFO(,,,,[])
AC_SUBST([LIBS])
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_LDFLAGS])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
