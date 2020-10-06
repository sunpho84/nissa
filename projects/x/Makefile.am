AM_CPPFLAGS=-I$(top_srcdir)/src

LDADD=$(top_builddir)/lib/libnissa.a

AUTOMAKE_OPTIONS=subdir-objects

x_sources= \
	src/diagrams/meson_exchange.cpp \
	src/diagrams/propagator_self_energy.cpp \
	src/diagrams/tadpole.cpp

x_sources+= \
	src/inverters/cg_Wilson_gluon_operator.cpp

x_sources+= \
	src/operators/Wilson_gluon_Klein_Gordon_operator.cpp \
	src/operators/twisted_Dirac_eoprec_operator.cpp

x_sources+= \
	src/propagators/twisted_propagator_g2_corr.cpp \
	src/propagators/Wilson_gluon_propagator.cpp

x_sources+= \
	src/routines/bmp.cpp \
	src/routines/correlations.cpp \
	src/routines/derivatives.cpp \
	src/routines/read_and_write.cpp \
	src/routines/shift.cpp

x_sources+= \
	src/stochastic/stochastic_source.cpp \
	src/stochastic/stochastic_tlSym_gluon_propagator.cpp \
	src/stochastic/stochastic_twisted_propagator.cpp 

x_sources+= \
	src/vertex/x_space_stochastic_qqg_vertex.cpp \
	src/vertex/vertex.cpp

##################################################### headers ###################################################

include_HEADERS= \
	src/diagrams/meson_exchange.hpp \
	src/diagrams/propagator_self_energy.hpp \
	src/diagrams/tadpole.hpp

include_HEADERS+= \
	src/inverters/cg_Wilson_gluon_operator.hpp

include_HEADERS+= \
	src/operators/Wilson_gluon_Klein_Gordon_operator.hpp \
	src/operators/twisted_Dirac_eoprec_operator.hpp

include_HEADERS+= \
	src/propagators/twisted_propagator_g2_corr.hpp \
	src/propagators/Wilson_gluon_propagator.hpp

include_HEADERS+= \
	src/routines/bmp.hpp \
	src/routines/correlations.hpp \
	src/routines/derivatives.hpp \
	src/routines/read_and_write.hpp \
	src/routines/shift.hpp

include_HEADERS+= \
	src/stochastic/stochastic_source.hpp \
	src/stochastic/stochastic_tlSym_gluon_propagator.hpp \
	src/stochastic/stochastic_twisted_propagator.hpp 

include_HEADERS+= \
	src/vertex/x_space_stochastic_qqg_vertex.hpp \
	src/vertex/vertex.hpp

############################################# define the executables ############################################

bin_PROGRAMS= \
	$(top_builddir)/test/x/twisted_propagator \
	$(top_builddir)/test/x/Wilson_gluon_propagator \
	$(top_builddir)/test/x/Wilson_gluon_Klein_Gordon_operator \
	$(top_builddir)/test/x/tlSym_gluon_stochastic_propagator \
	$(top_builddir)/test/x/tlSym_gluon_stochastic_propagator_no_null_mode \
	$(top_builddir)/test/x/stochastic_qqg_vertex \
	$(top_builddir)/test/x/twisted_propagator_stochastic_correction \
	$(top_builddir)/test/x/twisted_propagator_correction \
	$(top_builddir)/test/x/twisted_propagator_tadpole \
	$(top_builddir)/test/x/twisted_propagator_self \
	$(top_builddir)/test/x/stochastic_source_null_mode \
	$(top_builddir)/test/x/propagator_shift \
	$(top_builddir)/test/x/Fourier_transform \
	$(top_builddir)/test/x/meson_self_energy_stochastic_correction \
	$(top_builddir)/test/x/meson_exchange_correction \
	$(top_builddir)/test/x/meson_tad_pole_correction \
	$(top_builddir)/test/x/meson_tree_level_correction \
	$(top_builddir)/test/x/nazario \
	$(top_builddir)/test/x/nazario2

__top_builddir__test_x_twisted_propagator_SOURCES=$(x_sources) tests/test_twisted_propagator.cpp
__top_builddir__test_x_Wilson_gluon_propagator_SOURCES=$(x_sources) tests/test_Wilson_gluon_propagator.cpp
__top_builddir__test_x_Wilson_gluon_Klein_Gordon_operator_SOURCES=$(x_sources) tests/test_Wilson_gluon_Klein_Gordon_operator.cpp
__top_builddir__test_x_tlSym_gluon_stochastic_propagator_SOURCES=$(x_sources) tests/test_tlSym_gluon_stochastic_propagator.cpp
__top_builddir__test_x_tlSym_gluon_stochastic_propagator_no_null_mode_SOURCES=$(x_sources) tests/test_tlSym_gluon_stochastic_propagator_no_null_mode.cpp
__top_builddir__test_x_stochastic_qqg_vertex_SOURCES=$(x_sources) tests/test_stochastic_qqg_vertex.cpp
__top_builddir__test_x_twisted_propagator_stochastic_correction_SOURCES=$(x_sources) tests/test_twisted_propagator_stochastic_correction.cpp
__top_builddir__test_x_twisted_propagator_correction_SOURCES=$(x_sources) tests/test_twisted_propagator_correction.cpp
__top_builddir__test_x_twisted_propagator_tadpole_SOURCES=$(x_sources) tests/test_twisted_propagator_tadpole.cpp
__top_builddir__test_x_twisted_propagator_self_SOURCES=$(x_sources) tests/test_twisted_propagator_self.cpp
__top_builddir__test_x_stochastic_source_null_mode_SOURCES=$(x_sources) tests/test_stochastic_source_null_mode.cpp
__top_builddir__test_x_propagator_shift_SOURCES=$(x_sources) tests/test_propagator_shift.cpp
__top_builddir__test_x_Fourier_transform_SOURCES=$(x_sources) tests/test_Fourier_transform.cpp
__top_builddir__test_x_meson_self_energy_stochastic_correction_SOURCES=$(x_sources) tests/test_meson_self_energy_stochastic_correction.cpp
__top_builddir__test_x_meson_exchange_correction_SOURCES=$(x_sources) tests/test_meson_exchange_correction.cpp
__top_builddir__test_x_meson_tad_pole_correction_SOURCES=$(x_sources) tests/test_meson_tad_pole_correction.cpp
__top_builddir__test_x_meson_tree_level_correction_SOURCES=$(x_sources) tests/test_meson_tree_level_correction.cpp
__top_builddir__test_x_nazario_SOURCES=$(x_sources) tests/test_nazario.cpp
__top_builddir__test_x_nazario2_SOURCES=$(x_sources) tests/test_nazario2.cpp

bin_PROGRAMS+= \
	$(top_builddir)/bin/x/average_corr \
	$(top_builddir)/bin/x/bmp_scramble \
	$(top_builddir)/bin/x/compute_corrections \
	$(top_builddir)/bin/x/compute_tree_level_corrections \
	$(top_builddir)/bin/x/compute_self_energy_corrections \
	$(top_builddir)/bin/x/compute_tadpole_corrections \
	$(top_builddir)/bin/x/correct \
	$(top_builddir)/bin/x/combine \
	$(top_builddir)/bin/x/integrali_divergenti

__top_builddir__bin_x_average_corr_SOURCES=$(x_sources) progs/average_corr.cpp
__top_builddir__bin_x_bmp_scramble_SOURCES=$(x_sources) progs/bmp_scramble.cpp
__top_builddir__bin_x_combine_SOURCES=$(x_sources) progs/combine.cpp
__top_builddir__bin_x_compute_corrections_SOURCES=$(x_sources) progs/compute_corrections.cpp
__top_builddir__bin_x_compute_tree_level_corrections_SOURCES=$(x_sources) progs/compute_tree_level_corrections.cpp
__top_builddir__bin_x_compute_self_energy_corrections_SOURCES=$(x_sources) progs/compute_self_energy_corrections.cpp
__top_builddir__bin_x_compute_tadpole_corrections_SOURCES=$(x_sources) progs/compute_tadpole_corrections.cpp
__top_builddir__bin_x_integrali_divergenti_SOURCES=$(x_sources) progs/integrali_divergenti.cpp
__top_builddir__bin_x_correct_SOURCES=$(x_sources) progs/correct.cpp
