include $(top_srcdir)/config/baseMakefile.am

AM_CPPFLAGS=-I$(top_srcdir)/src

AUTOMAKE_OPTIONS=subdir-objects

bin_PROGRAMS+= \
	$(top_builddir)/tools/conf_measure/topo \
	$(top_builddir)/tools/conf_transform/convert/eo_to_ildg \
	$(top_builddir)/tools/conf_transform/convert/gpu_to_ildg \
	$(top_builddir)/tools/conf_transform/convert/ildg_to_eo \
        $(top_builddir)/tools/conf_transform/convert/ildg_to_Claudio \
        $(top_builddir)/tools/conf_transform/convert/Claudio_to_ildg \
	$(top_builddir)/tools/conf_transform/convert/ildg_to_gpu \
	$(top_builddir)/tools/conf_transform/create/create_cold \
	$(top_builddir)/tools/conf_transform/create/create_hot \
	$(top_builddir)/tools/conf_transform/gaugefix/landau \
	$(top_builddir)/tools/conf_transform/gaugefix/random \
	$(top_builddir)/tools/conf_transform/smear/hex \
	$(top_builddir)/tools/conf_transform/smear/stout \
	$(top_builddir)/tools/conf_transform/edit/build_from_raw \
	$(top_builddir)/tools/conf_transform/edit/grab_ildg_data \
	$(top_builddir)/tools/conf_transform/edit/rewind_random \
	$(top_builddir)/tools/conf_transform/edit/skip_random \
	$(top_builddir)/tools/conf_transform/edit/unitarize \
	$(top_builddir)/tools/find_best_vir_partitioning/find \
	$(top_builddir)/tools/endianness/endianness_check \
	$(top_builddir)/tools/endianness/endianness_swap \
	$(top_builddir)/tools/Gaussiansmear_test \
	$(top_builddir)/tools/print_gamma/gamma_test \
	$(top_builddir)/tools/remez_algorithm/remez_create \
	$(top_builddir)/tools/read_ildg/read_ildg \
	$(top_builddir)/tools/read_ildg/read_double_vector \
	$(top_builddir)/tools/spinors/ildg_to_txt \
	$(top_builddir)/tools/spinors/take_parity \
	$(top_builddir)/tools/spinors/sources_compare \
	$(top_builddir)/tools/unitarity_check/unitarity_check

#__top_builddir__tools_benchmarks_stag_SOURCES=benchmarks/stag.cpp
__top_builddir__tools_find_best_vir_partitioning_find_SOURCES=find_best_vir_partitioning/find.cpp
__top_builddir__tools_remez_algorithm_remez_create_SOURCES=remez_algorithm/remez_create.cpp
__top_builddir__tools_conf_measure_topo_SOURCES=conf_measure/topo.cpp
__top_builddir__tools_conf_transform_convert_gpu_to_ildg_SOURCES=conf_transform/convert/gpu_to_ildg.cpp conf_transform/convert/gpu_stagphase.cpp
__top_builddir__tools_conf_transform_convert_ildg_to_eo_SOURCES=conf_transform/convert/ildg_to_eo.cpp
__top_builddir__tools_conf_transform_convert_Claudio_to_ildg_SOURCES=conf_transform/convert/Claudio_to_ildg.cpp
__top_builddir__tools_conf_transform_convert_ildg_to_Claudio_SOURCES=conf_transform/convert/ildg_to_Claudio.cpp
__top_builddir__tools_conf_transform_convert_eo_to_ildg_SOURCES=conf_transform/convert/eo_to_ildg.cpp
__top_builddir__tools_conf_transform_convert_ildg_to_gpu_SOURCES=conf_transform/convert/ildg_to_gpu.cpp conf_transform/convert/gpu_stagphase.cpp
__top_builddir__tools_conf_transform_create_create_cold_SOURCES=conf_transform/create/create_cold.cpp
__top_builddir__tools_endianness_endianness_check_SOURCES=endianness/endianness_check.cpp
__top_builddir__tools_endianness_endianness_swap_SOURCES=endianness/endianness_swap.cpp
__top_builddir__tools_read_ildg_read_ildg_SOURCES=read_ildg/read_ildg.cpp
__top_builddir__tools_read_ildg_read_double_vector_SOURCES=read_ildg/read_double_vector.cpp
__top_builddir__tools_conf_transform_create_create_hot_SOURCES=conf_transform/create/create_hot.cpp
__top_builddir__tools_conf_transform_gaugefix_landau_SOURCES=conf_transform/gaugefix/landau.cpp
__top_builddir__tools_conf_transform_gaugefix_random_SOURCES=conf_transform/gaugefix/random.cpp
__top_builddir__tools_conf_transform_smear_hex_SOURCES=conf_transform/smear/hex.cpp
__top_builddir__tools_conf_transform_smear_stout_SOURCES=conf_transform/smear/stout.cpp
__top_builddir__tools_conf_transform_edit_build_from_raw_SOURCES=conf_transform/edit/build_from_raw.cpp
__top_builddir__tools_conf_transform_edit_grab_ildg_data_SOURCES=conf_transform/edit/grab_ildg_data.cpp
__top_builddir__tools_conf_transform_edit_rewind_random_SOURCES=conf_transform/edit/rewind_random.cpp
__top_builddir__tools_conf_transform_edit_skip_random_SOURCES=conf_transform/edit/skip_random.cpp
__top_builddir__tools_conf_transform_edit_unitarize_SOURCES=conf_transform/edit/unitarize.cpp
__top_builddir__tools_Gaussiansmear_test_SOURCES=Gaussiansmear_test.cpp
__top_builddir__tools_spinors_ildg_to_txt_SOURCES=spinors/ildg_to_txt.cpp
__top_builddir__tools_spinors_sources_compare_SOURCES=spinors/sources_compare.cpp
__top_builddir__tools_spinors_take_parity_SOURCES=spinors/take_parity.cpp
__top_builddir__tools_print_gamma_gamma_test_SOURCES=print_gamma/gamma_test.cpp
__top_builddir__tools_unitarity_check_unitarity_check_SOURCES=unitarity_check/unitarity_check.cpp

EXTRA_DIST+=find_best_vir_partitioning/prepare_mapping.sh

#generate assembly report if asked
if ASSEMBLY_REPORT
BUILT_SOURCES+=$(assembly_reports)
CLEANFILES+=$(assembly_reports)
endif

#	$(top_builddir)/tools/benchmarks/stag
