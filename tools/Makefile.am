AM_CPPFLAGS=-I$(top_srcdir)/src
LDADD=$(top_builddir)/lib/libnissa.la
AUTOMAKE_OPTIONS=subdir-objects

bin_PROGRAMS= \
	$(top_builddir)/tools/find_best_vir_partitioning/find \
	$(top_builddir)/tools/conf_transform/gpu_to_ildg \
	$(top_builddir)/tools/conf_transofrm/ildg_to_gpu \
	$(top_builddir)/tools/conf_transform/smear

__top_builddir__tools_find_best_vir_partitioning_find_SOURCES=find_best_vir_partitioning/find.cpp
__top_builddir__tools_conf_transform_gpu_to_ildg_SOURCES=conf_transform/gpu_to_ildg.cpp
__top_builddir__tools_conf_transform_ildg_to_gpu_SOURCES=conf_transform/ildg_to_gpu.cpp
__top_builddir__tools_conf_transform_smear_SOURCES=conf_transform/smear.cpp