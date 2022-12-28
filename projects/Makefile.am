AUTOMAKE_OPTIONS=subdir-objects

AM_LFLAGS = -s -P`basename $* | sed 's,l$$,,'`_ -olex.yy.c
AM_YFLAGS = -Wno-yacc -d -p `basename $* | sed 's,y$$,,'`_

LDADD=$(top_builddir)/lib/libnissa.a
AM_CPPFLAGS=-I$(top_srcdir)/src

#init to empty
bin_PROGRAMS=
include_HEADERS=
BUILT_SOURCES=
CLEANFILES=

#include $(srcdir)/generate_confs/Makefile.am
include $(srcdir)/ib/Makefile.am
#include $(srcdir)/nucleons/Makefile.am
#include $(srcdir)/reno_const/Makefile.am
#include $(srcdir)/semileptonic/Makefile.am

if NDIM2
include $(srcdir)/cpn/Makefile.am
endif
