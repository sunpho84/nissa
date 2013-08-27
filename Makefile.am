#include m4 dir
ACLOCAL_AMFLAGS=-I m4

#include subdirs
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
	src/routines \
	projects

#add nissa library
lib_LTLIBRARIES=lib/libnissa.la
lib_libnissa_la_SOURCES=
lib_libnissa_la_LIBADD= \
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

#adds extra files to distribution
EXTRA_DIST= \
        svn_revision \
        README \
        src/dirac_operators/dirac_operator_WclovQ/dirac_operator_WclovQ_portable.cpp \
        src/dirac_operators/dirac_operator_stD/dirac_operator_stD_portable.cpp \
        src/dirac_operators/dirac_operator_tmDeoimpr/dirac_operator_tmDeoimpr_portable.cpp \
        src/dirac_operators/dirac_operator_tmQ/dirac_operator_tmQ_portable.cpp \
        src/dirac_operators/dirac_operator_tmQ_left/dirac_operator_tmQ_left_portable.cpp \
        src/hmc/rootst_eoimpr/rat_expansion_database.cpp \
        src/inverters/twisted_mass/cg_invert_tmQ2_common.cpp \
        src/nissa.h \
        src/create_headers.sh \
        src/inverters/templates/cg_128_invert_template.cpp \
        src/inverters/templates/cg_invert_template.cpp \
        src/inverters/templates/cgm_invert_template.cpp \
        src/inverters/templates/cg_128_invert_template_threaded.cpp \
        src/inverters/templates/cg_invert_template_threaded.cpp \
        src/inverters/templates/cgm_invert_template_threaded.cpp \
        tools/clusterize.cpp \
        tools/clusterize2.cpp \
        tools/pack.c
