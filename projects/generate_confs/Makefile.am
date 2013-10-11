AM_CPPFLAGS=-I$(top_srcdir)/src

include $(top_srcdir)/config/cuda.mk

AUTOMAKE_OPTIONS=subdir-objects

LDADD=$(top_builddir)/lib/libnissa.la $(CUDA_LIBS)

bin_PROGRAMS=$(top_builddir)/bin/generate_stag 

__top_builddir__bin_generate_stag_SOURCES=generate_stag.cpp
__top_builddir__bin_generate_stag_LDADD=$(LDADD)
__top_builddir__bin_generate_stag_LDFLAGS=$(LDFLAGS) $(CUDA_FLAGS)
