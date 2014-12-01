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
	src/create_headers.sh \
	tools/ave.cpp \
	tools/common_tools.cpp \
	tools/clusterize.cpp \
	tools/clusterize2.cpp \
	tools/pack.c

DIST_SUBDIRS=${SUBDIRS} tools
