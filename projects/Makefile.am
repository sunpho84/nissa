AUTOMAKE_OPTIONS=subdir-objects

AM_LFLAGS = -s -P`basename $* | sed 's,l$$,,'`_ -olex.yy.c
AM_YFLAGS = -Wno-yacc -d -p `basename $* | sed 's,y$$,,'`_

AM_CPPFLAGS=-I$(top_srcdir)/src

include $(top_srcdir)/config/baseMakefile.am

include $(srcdir)/ib/Makefile.am
