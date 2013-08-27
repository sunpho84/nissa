ACLOCAL_AMFLAGS=-I m4
SUBDIRS= \
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
	src/routines
lib_LTLIBRARIES=libnissa.la
libnissa_la_SOURCES=
libnissa_la_LIBADD= \
	src/base/libbase.la \
	src/communicate/libcommunicate.la \
	src/dirac_operators/libdirac_operators.la \
	src/geometry/libgeometry.la \
	src/hmc/libhmc.la \
	src/inverters/libinverters.la \
	src/io/libio.la \
	src/new_types/libnew_types.la \
	src/operations/liboperations.la \
	src/routines/libroutines.la