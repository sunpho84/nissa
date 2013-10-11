include_HEADERS=nissa.hpp

#add nissa library
lib_LTLIBRARIES=$(top_builddir)/lib/libnissa.la
__top_builddir__lib_libnissa_la_SOURCES=
__top_builddir__lib_libnissa_la_LIBADD= \
        base/libbase.la \
        communicate/libcommunicate.la \
        dirac_operators/libdirac_operators.la \
        geometry/libgeometry.la \
        hmc/libhmc.la \
        inverters/libinverters.la \
        io/libio.la \
        linalgs/liblinalgs.la \
        new_types/libnew_types.la \
        operations/liboperations.la \
        routines/libroutines.la

#include subdirs
SUBDIRS= \
        base \
        communicate \
        dirac_operators \
        geometry \
        hmc \
        inverters \
        io \
        linalgs \
        new_types \
        operations \
        routines

if GPU
SUBDIRS+=gpu
__top_builddir__lib_libnissa_la_LIBADD+=gpu/libgpu.la
endif

if BGQ
SUBDIRS+=bgq
__top_builddir__lib_libnissa_la_LIBADD+=bgq/libbgq.la
endif

