# $Id$
## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS = -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\" 

lib_LTLIBRARIES = libgridcount@LIBSUFFIX@.la
libgridcount@LIBSUFFIX@_la_LDFLAGS        = -version-info @SHARED_VERSION_INFO@


libgridcount@LIBSUFFIX@_la_SOURCES = \
        utilgmx.c utilgmx.h   xf.c xf.h   plt.c plt.h \
        count.c count.h   grid3D.c grid3D.h   xdr_grid.c xdr_grid.h \
	gridcount.h


bin_PROGRAMS = 	g_ri3Dc a_ri3Dc

g_ri3Dc_SOURCES	= g_ri3Dc.c
a_ri3Dc_SOURCES	= a_ri3Dc.c  


LDADD = $(lib_LTLIBRARIES) ./libgridcount@LIBSUFFIX@.la \
	../../mdlib/libmd@LIBSUFFIX@.la \
	../../gmxlib/libgmx@LIBSUFFIX@.la 


# link the mpi library to non-mpi names if the latter are not present
install-exec-hook:
	libname="libgridcount@LIBSUFFIX@"; \
	nompi="`echo $$libname | sed -e 's,_mpi,,'`"; \
	libdir="$(libdir)"; \
	if echo $$libname | grep mpi >/dev/null ; then \
	  (cd $$libdir && test -e $$libname.a -a ! -e $$nompi.a && $(LN_S) $$libname.a $$nompi.a ; exit 0); \
	  (cd $$libdir && test -e $$libname.so -a ! -e $$nompi.so && $(LN_S) $$libname.so $$nompi.so ; exit 0); \
	fi; 

CLEANFILES   = *.la *~ \\\#*














