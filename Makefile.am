include $(top_srcdir)/config/baseMakefile.am

#include m4 dir
ACLOCAL_AMFLAGS=-I m4

#include subdirs
SUBDIRS= \
	src \
	. \
	projects

#adds extra files to distribution
EXTRA_DIST= \
	config \
	README \
	src/nissa.hpp \
	tools/ave.cpp \
	tools/common_tools.cpp \
        tools/conf_transform \
	tools/clusterize.cpp \
	tools/clusterize2.cpp \
	tools/pack.c

DIST_SUBDIRS=${SUBDIRS} tools

#generate assembly report if asked
if ASSEMBLY_REPORT
BUILT_SOURCES+=$(assembly_reports)
CLEANFILES+=$(assembly_reports)
endif
