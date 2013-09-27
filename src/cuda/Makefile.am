AM_CPPFLAGS=-I$(top_srcdir)/src

include $(top_srcdir)/config/cuda.mk

LINK=$(LIBTOOL) --mode=link $(CC) -o $@ $(CUDA_LDFLAGS) $(CUDA_LIBS)

noinst_LTLIBRARIES=libcuda.la
libcuda_la_SOURCES= \
	fields.cu \
	global_variables.cu \
	init.cu \
	wrappers.cu

cudaincludedir=$(includedir)/cuda
cudainclude_HEADERS= \
	fields.hpp \
	global_variables.hpp \
	init.hpp \
	macros.hpp \
	new_types.hpp \
        wrappers.hpp
