AM_CPPFLAGS=-I$(top_srcdir)/src

#library name
noinst_LTLIBRARIES=libbase.la

#sources
libbase_la_SOURCES= \
	bench.cpp \
	close.cpp \
	debug.cpp \
	global_variables.cpp \
	init.cpp \
	random.cpp \
	vectors.cpp

#headers
baseincludedir=$(includedir)/base
baseinclude_HEADERS= \
	bench.hpp \
	close.hpp \
	debug.hpp \
	global_variables.hpp \
	init.hpp \
	macros.hpp \
	random.hpp \
	thread_macros.hpp \
	vectors.hpp

#include bgq files
if BGQ
baseinclude_HEADERS+=bgq_instructions.hpp
endif

#svn versioning: specify to build sources and include script
BUILT_SOURCES=update_svnversion svnversion.hpp
EXTRA_DIST=update_svnversion.hpp.sh
update_svnversion:
	$(srcdir)/update_svnversion.hpp.sh $(top_srcdir)
.PHONY: update_svnversion

#include global_variables.cpp in the installation
#needed because included in global_variables.h
global_variabledir=$(prefix)/include/base
global_variable_DATA=global_variables.cpp

CLEANFILES=svnversion.hpp