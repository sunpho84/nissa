AM_CXXFLAGS=-I$(top_srcdir)/src

BUILT_SOURCES=update_svnversion svnversion.h

if BGQ
base_bgq_header_files=../bgq/spi.h bgq_instructions.h
endif

noinst_LTLIBRARIES=libbase.la
libbase_la_SOURCES= \
	bench.cpp \
	close.cpp \
	debug.cpp \
	global_variables.cpp \
	init.cpp \
	random.cpp \
	vectors.cpp

baseincludedir=$(includedir)/base
baseinclude_HEADERS= \
	$(base_bgq_header_files) \
	bench.h \
	close.h \
	debug.h \
	global_variables.h \
	init.h \
	macros.h \
	random.h \
	thread_macros.h \
	vectors.h

EXTRA_DIST=update_svnversion.h.sh svnversion.h

update_svnversion:
	cd $(srcdir) && ./update_svnversion.h.sh

.PHONY: update_svnversion