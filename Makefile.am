#include m4 dir
ACLOCAL_AMFLAGS=-I m4

if BGQ
bgq_subdirs=src/bgq
bgq_sublib=src/bgq/libbgq.la
else
bgq_subdirs=
bgq_sublib=
endif

#include subdirs
SUBDIRS= \
	$(bgq_subdirs) \
	src/base \
	src/communicate \
	src/dirac_operators \
	src/geometry \
	src/hmc \
	src/inverters \
	src/io \
	src/linalgs \
	src/new_types \
	src/operations \
	src/routines \
	. \
	projects

#add nissa library
lib_LTLIBRARIES=lib/libnissa.la
lib_libnissa_la_SOURCES=
lib_libnissa_la_LIBADD= \
	$(bgq_sublib) \
	src/base/libbase.la \
	src/communicate/libcommunicate.la \
	src/dirac_operators/libdirac_operators.la \
	src/geometry/libgeometry.la \
	src/hmc/libhmc.la \
	src/inverters/libinverters.la \
	src/io/libio.la \
	src/linalgs/liblinalgs.la \
	src/new_types/libnew_types.la \
	src/operations/liboperations.la \
	src/routines/libroutines.la

#adds extra files to distribution
EXTRA_DIST= \
	libtool \
	config \
	bootstrap \
        README \
        src/nissa.h \
        src/create_headers.sh \
        tools/clusterize.cpp \
        tools/clusterize2.cpp \
        tools/pack.c
