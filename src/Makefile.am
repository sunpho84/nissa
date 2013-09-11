if BGQ
bgq_subdirs=bgq
bgq_sublib=bgq/libbgq.la
else
bgq_subdirs=
bgq_sublib=
endif

include_HEADERS=nissa.h

#add nissa library
lib_LTLIBRARIES=$(top_builddir)/lib/libnissa.la
__top_builddir__lib_libnissa_la_SOURCES=
__top_builddir__lib_libnissa_la_LIBADD= \
        base/libbase.la \
        $(bgq_sublib) \
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
        $(bgq_subdirs) \
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


