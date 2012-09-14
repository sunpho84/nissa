#externally defined path and variables
CC=mpicxx
CFLAGS=-O5

#derived parameters
SVN_VERS:=$(shell svnversion 2>/dev/null)
MACROS= #$(addprefix -D,SVN_VERS=\"$(SVN_VERS)\")
INCLUDE_PATH=../../src
LIBRARY_PATH=

GCC=gcc

################################################ define the objects ############################################

#include all objects
objects=$(addprefix src/, \
	$(addprefix diagrams/, meson_exchange propagator_self_energy tadpole) \
	$(addprefix inverters/, cg_eoprec_twisted_free_operator cg_Wilson_gluon_operator) \
	$(addprefix operators/, Wilson_gluon_Klein_Gordon_operator twisted_Dirac_eoprec_operator) \
	$(addprefix propagators/, twisted_propagator twisted_propagator_g2_corr Wilson_gluon_propagator tlSym_gluon_propagator) \
	$(addprefix routines/, bmp correlations derivatives fourier read_and_write shift) \
	$(addprefix stochastic/, stochastic_source stochastic_tlSym_gluon_propagator stochastic_twisted_propagator) \
	$(addprefix types/, types_routines) \
	$(addprefix vertex/, x_space_stochastic_qqg_vertex vertex))

############################################# define the executables ############################################

#test programs
tests=$(addprefix tests/, \
	test_twisted_propagator \
	test_Wilson_gluon_propagator \
	test_tlSym_gluon_propagator \
	test_Wilson_gluon_Klein_Gordon_operator \
	test_tlSym_gluon_stochastic_propagator \
	test_tlSym_gluon_stochastic_propagator_no_null_mode \
	test_stochastic_qqg_vertex \
	test_twisted_propagator_stochastic_correction \
	test_twisted_propagator_correction \
	test_stochastic_source_null_mode \
	test_propagator_shift \
	test_Fourier_transform \
	test_meson_self_energy_stochastic_correction \
	test_meson_exchange_correction \
	test_meson_tree_level_correction \
	test_nazario \
	test_nazario2 \
)

#full programs
progs=$(addprefix progs/, \
	compute_corrections \
	compute_tree_level_corrections \
	compute_self_energy_corrections \
	compute_tadpole_corrections \
	average_corr \
	bmp_scramble \
	correct \
	combine \
)

#collect all the executables
executables=$(tests) $(progs)

################################################## global targets ###############################################

all: Makefile $(executables)

clean:
	 rm -rf $(addsuffix .d,$(objects) $(executables)) $(addsuffix .o,$(objects)) $(executables)

#include table of dependencies
-include $(addsuffix .d,$(projects) $(executables))

############################################ rules to produce objects ############################################

$(addsuffix .o,$(objects) $(executables)): %.o: %.cpp %.d Makefile
	$(CC)				\
	$(addprefix -I,$(INCLUDE_PATH)) \
	$<				\
	$(CFLAGS) 			\
	$(MACROS) 			\
	-c -o $@

######################################### rules to produce dependencies ##########################################

$(addsuffix .d,$(objects) $(executables)): %.d: %.cpp Makefile
	$(GCC) $< -MM -MT $(@:%.d=%.o) $(addprefix -I,$(INCLUDE_PATH)) -o $@

######################################### rules to produce executables ############################################

$(executables): %: %.o $(addsuffix .o,$(objects)) ../../libnissa.a Makefile
	$(CC)                           \
	$(addprefix -I,$(INCLUDE_PATH)) \
	$(addprefix -L,$(LIBRARY_PATH)) \
	$(CFLAGS)                       \
	$(MACROS)                       \
	-o $@                           \
	$(addsuffix .o, $@)           \
	$(addsuffix .o, $(objects))     \
	../../libnissa.a

############################################## phonyfyse the needed target ##########################################

.PHONY: clean
