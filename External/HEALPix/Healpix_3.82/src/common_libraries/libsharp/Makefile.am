ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libsharp.la

libsharp_la_SOURCES = \
  c_utils/c_utils.c \
  c_utils/c_utils.h \
  c_utils/walltime_c.c \
  c_utils/walltime_c.h \
  pocketfft/pocketfft.c \
  pocketfft/pocketfft.h \
  libsharp/sharp.c \
  libsharp/sharp_almhelpers.c \
  libsharp/sharp_core.c \
  libsharp/sharp_geomhelpers.c \
  libsharp/sharp_legendre_roots.c \
  libsharp/sharp_ylmgen_c.c \
  libsharp/sharp_internal.h \
  libsharp/sharp_legendre_roots.h \
  libsharp/sharp_vecsupport.h \
  libsharp/sharp_ylmgen_c.h

if HAVE_MULTIARCH

libavx_la_SOURCES = libsharp/sharp_core_inc.c
libavx2_la_SOURCES = libsharp/sharp_core_inc.c
libfma_la_SOURCES = libsharp/sharp_core_inc.c
libfma4_la_SOURCES = libsharp/sharp_core_inc.c
libavx512f_la_SOURCES = libsharp/sharp_core_inc.c

noinst_LTLIBRARIES = libavx.la libavx2.la libfma.la libfma4.la libavx512f.la

libavx_la_CFLAGS = ${AM_CFLAGS} -mavx -DARCH=avx
libavx2_la_CFLAGS = ${AM_CFLAGS} -mavx2 -DARCH=avx2
libfma_la_CFLAGS = ${AM_CFLAGS} -mfma -DARCH=fma
libfma4_la_CFLAGS = ${AM_CFLAGS} -mfma4 -DARCH=fma4
libavx512f_la_CFLAGS = ${AM_CFLAGS} -mavx512f -DARCH=avx512f

libsharp_la_LIBADD = libavx.la libavx2.la libfma.la libfma4.la libavx512f.la

endif

# format is "current:revision:age"
# - Always increase the revision value.
# - Increase the current value whenever an interface has been added, removed or changed.
# - Increase the age value only if the changes made to the ABI are backward compatible.
libsharp_la_LDFLAGS = -version-info 2:2:0

nobase_include_HEADERS = \
  libsharp/sharp.h \
  libsharp/sharp_geomhelpers.h \
  libsharp/sharp_almhelpers.h \
  libsharp/sharp_cxx.h

EXTRA_DIST = \
  runtest.sh \
  COMPILE

check_PROGRAMS = sharp_testsuite
sharp_testsuite_SOURCES = libsharp/sharp_testsuite.c c_utils/memusage.c c_utils/memusage.h
sharp_testsuite_LDADD = libsharp.la

TESTS = runtest.sh

AM_CFLAGS = @AM_CFLAGS@

pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = @PACKAGE_NAME@.pc

DISTCLEANFILES=@PACKAGE_NAME@.pc @PACKAGE_NAME@.pc.in  @PACKAGE_NAME@-uninstalled.pc  @PACKAGE_NAME@-uninstalled.sh
