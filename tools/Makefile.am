AM_CPPFLAGS=-I$(top_srcdir)/src
LDADD=$(top_builddir)/lib/libnissa.a
AUTOMAKE_OPTIONS=subdir-objects

bin_PROGRAMS= \
	$(top_builddir)/tools/find_best_vir_partitioning/find \
	$(top_builddir)/tools/conf_measure/topo \
	$(top_builddir)/tools/conf_transform/convert/eo_to_ildg \
	$(top_builddir)/tools/conf_transform/convert/gpu_to_ildg \
	$(top_builddir)/tools/conf_transform/convert/ildg_to_eo \
	$(top_builddir)/tools/conf_transform/convert/ildg_to_gpu \
	$(top_builddir)/tools/conf_transform/create/create_cold \
	$(top_builddir)/tools/conf_transform/create/create_hot \
	$(top_builddir)/tools/conf_transform/smear/stout \
	$(top_builddir)/tools/conf_transform/edit/grab_ildg_data \
	$(top_builddir)/tools/conf_transform/edit/skip_random \
	$(top_builddir)/tools/print_gamma/gamma_test \
	$(top_builddir)/tools/remez_algorithm/remez_create \
	$(top_builddir)/tools/unitarity_check/unitarity_check \
	$(top_builddir)/tools/x_clusterize/x_clusterize

__top_builddir__tools_find_best_vir_partitioning_find_SOURCES=find_best_vir_partitioning/find.cpp
__top_builddir__tools_remez_algorithm_remez_create_SOURCES=remez_algorithm/remez_create.cpp
__top_builddir__tools_conf_measure_topo_SOURCES=conf_measure/topo.cpp
__top_builddir__tools_conf_transform_convert_gpu_to_ildg_SOURCES=conf_transform/convert/gpu_to_ildg.cpp
__top_builddir__tools_conf_transform_convert_ildg_to_eo_SOURCES=conf_transform/convert/ildg_to_eo.cpp
__top_builddir__tools_conf_transform_convert_eo_to_ildg_SOURCES=conf_transform/convert/eo_to_ildg.cpp
__top_builddir__tools_conf_transform_convert_ildg_to_gpu_SOURCES=conf_transform/convert/ildg_to_gpu.cpp
__top_builddir__tools_conf_transform_create_create_cold_SOURCES=conf_transform/create/create_cold.cpp
__top_builddir__tools_conf_transform_create_create_hot_SOURCES=conf_transform/create/create_hot.cpp
__top_builddir__tools_conf_transform_smear_stout_SOURCES=conf_transform/smear/stout.cpp
__top_builddir__tools_conf_transform_edit_grab_ildg_data_SOURCES=conf_transform/edit/grab_ildg_data.cpp
__top_builddir__tools_conf_transform_edit_skip_random_SOURCES=conf_transform/edit/skip_random.cpp
__top_builddir__tools_print_gamma_gamma_test_SOURCES=print_gamma/gamma_test.cpp
__top_builddir__tools_unitarity_check_unitarity_check_SOURCES=unitarity_check/unitarity_check.cpp
__top_builddir__tools_x_clusterize_x_clusterize_SOURCES=x_clusterize/x_clusterize.cpp

EXTRA_DIST=find_best_vir_partitioning/prepare_mapping.sh
