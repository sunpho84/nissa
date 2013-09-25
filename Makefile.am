#include m4 dir
ACLOCAL_AMFLAGS=-I m4

#include cude
if CUDA
include $(top_srcdir)/config/cuda.mk
endif

#include subdirs
SUBDIRS= \
	src \
	. \
	projects

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

#remove libtool from distribution
dist-hook:
	rm -fr $(top_distdir)/libtool