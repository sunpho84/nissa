AM_CPPFLAGS=-I$(top_srcdir)/src
LDADD=$(top_builddir)/lib/libnissa.la
AUTOMAKE_OPTIONS=subdir-objects

bin_PROGRAMS= \
	find_best_vir_partitioning/find \
	conf_convert/gpu_to_ildg \
	conf_convert/ildg_to_gpu

find_best_vir_partitioning_find_SOURCES=find_best_vir_partitioning/find.cpp
conf_convert_gpu_to_ildg_SOURCES=conf_convert/gpu_to_ildg.cpp
conf_convert_ildg_to_gpu_SOURCES=conf_convert/ildg_to_gpu.cpp