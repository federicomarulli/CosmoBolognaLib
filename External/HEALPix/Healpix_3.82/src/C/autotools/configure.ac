AC_INIT([chealpix], [3.60.0])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

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

AC_PROG_CC_C99

dnl
dnl Determine the compiler flags (CFITSIO_CFLAGS)
dnl and linker flags (CFITSIO_LIBS) required to build against libcfitsio.
dnl
dnl Add help messages to describe the variables CFITSIO_CFLAGS and CFITSIO_LIBS,
dnl which can be provided as environment variables or as command-line arguments.
dnl These are output variables of PKG_CHECK_MODULES, but they can be overridden.
dnl
AC_ARG_VAR([CFITSIO_CFLAGS], [compiler flags needed for building with cfitsio, if not found by pkg-config])
AC_ARG_VAR([CFITSIO_LIBS], [linker flags needed for building with cfitsio, if not found by pkg-config])
dnl
dnl Now check for cfitsio using pkg-config. If cfitsio is not found with
dnl pkg-config, CFITSIO_LIBS defaults to "-lcfitsio".
dnl
PKG_CHECK_MODULES([CFITSIO], [cfitsio], , [CFITSIO_LIBS=-lcfitsio])
dnl
dnl Linking against libcfitsio requires libm too.
dnl
CFITSIO_LIBS="$CFITSIO_LIBS -lm"
dnl
dnl With pkg-config < 0.24, PKG_CHECK_MODULES does not AC_SUBST the _CFLAGS
dnl and _LIBS variables.
dnl
AC_SUBST([CFITSIO_CFLAGS])
AC_SUBST([CFITSIO_LIBS])
dnl
dnl Perform a test compile and link to determine if cfitsio is usable. If
dnl pkg-config did not find cfitsio, then the user may still have provided
dnl CFITSIO_CFLAGS and CFITSIO_LIBS manually. If pkg-config did find cfitsio,
dnl it does not hurt to test that cfitsio is usable.
dnl
save_CFLAGS=$CFLAGS
save_LIBS=$LIBS
CFLAGS="$CFLAGS $CFITSIO_CFLAGS"
LIBS="$LIBS $CFITSIO_LIBS"
AC_CHECK_LIB([cfitsio], [ffgnrwll], [:], AC_MSG_ERROR([could not find the cfitsio library]))
AC_CHECK_HEADERS([fitsio.h], [:], AC_MSG_ERROR([could not find the cfitsio header file]))
CFLAGS=$save_CFLAGS
LIBS=$save_LIBS

dnl This macro affects conditional compilation inside some of the source files.
dnl FIXME: do we want to make cfitsio optional? It is required for healpix_cxx.
AC_DEFINE([ENABLE_FITSIO])

AC_PROG_LIBTOOL

dnl
dnl Create pkgconfig .pc file.
dnl
AX_CREATE_PKGCONFIG_INFO(,,,,[$CFITSIO_CFLAGS])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
