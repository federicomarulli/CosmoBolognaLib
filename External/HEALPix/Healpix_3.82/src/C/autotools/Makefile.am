ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libchealpix.la

include_HEADERS = chealpix.h

libchealpix_la_SOURCES = chealpix.c chealpix.h
libchealpix_la_LIBADD = $(CFITSIO_LIBS)

TESTS = test_chealpix
check_PROGRAMS = $(TESTS)
test_chealpix_LDADD = libchealpix.la

AM_CPPFLAGS = $(CFITSIO_CFLAGS)

pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = @PACKAGE_NAME@.pc

DISTCLEANFILES=@PACKAGE_NAME@.pc @PACKAGE_NAME@.pc.in  @PACKAGE_NAME@-uninstalled.pc  @PACKAGE_NAME@-uninstalled.sh
CLEANFILES = test_output.fits
