MAINTAINERCLEANFILES = Makefile.in
AM_CPPFLAGS = \
-I. \
-I$(top_srcdir)/src/lib \
@GSL_CFLAGS@ \
@BLAS_CFLAGS@ \
@POPT_CFLAGS@ \
@IPP_CFLAGS@ 

libSPH = $(top_builddir)/src/lib/.libs/libSPH.la \
@GSL_LIBS@ \
@BLAS_LIBS@ \
@POPT_LIBS@ \
@IPP_LIBS@ \
-lpthread

bin_PROGRAMS = sph

sph_LDADD = $(libSPH)
sph_SOURCES = sph.c

