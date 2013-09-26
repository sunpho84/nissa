AM_CPPFLAGS=-I$(top_srcdir)/src

BUILT_SOURCES=update_svnversion svnversion.hpp

if BGQ
base_bgq_header_files=../bgq/spi.hpp bgq_instructions.hpp
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
	bench.hpp \
	close.hpp \
	debug.hpp \
	global_variables.hpp \
	init.hpp \
	macros.hpp \
	random.hpp \
	thread_macros.hpp \
	vectors.hpp

EXTRA_DIST=update_svnversion.hpp.sh svnversion.hpp

#include global_variables.cpp in the installation
global_variabledir=$(prefix)/include/base
global_variable_DATA=global_variables.cpp

update_svnversion:
	cd $(srcdir) && ./update_svnversion.hpp.sh

.PHONY: update_svnversion
