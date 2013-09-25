AUTOMAKE_OPTIONS=subdir-objects
include $(top_srcdir)/config/cuda.mk
AM_CXXFLAGS=-I$(top_srcdir)/src

LDADD=$(top_builddir)/lib/libnissa.la $(CUDA_LIBS)

bin_PROGRAMS=$(top_builddir)/bin/generate_stag

__top_builddir__bin_generate_stag_SOURCES=generate_stag.cpp
__top_builddir__bin_generate_stag_LDADD=$(LDADD)
__top_builddir__bin_generate_stag_LDFLAGS=$(LDFLAGS) $(CUDA_LDFLAGS)
