AM_CPPFLAGS=-I$(top_srcdir)/src

include $(top_srcdir)/config/cuda.mk

LINK=$(LIBTOOL) --mode=link $(CC) -o $@ $(CUDA_LDFLAGS) $(CUDA_LIBS)

noinst_LTLIBRARIES=libgpu.la
libgpu_la_SOURCES= \
	fields.cu \
	global_variables.cu \
	init.cu \
	new_types.cu \
	wrappers.cu

gpuincludedir=$(includedir)/gpu
gpuinclude_HEADERS= \
	fields.hpp \
	global_variables.hpp \
	init.hpp \
	macros.hpp \
	new_types.hpp \
        wrappers.hpp
