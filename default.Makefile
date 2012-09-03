#externally defined path and variables
CC=SED_CC
CFLAGS=SED_CFLAGS

#derived parameters
SVN_VERS:=$(shell svnversion 2>/dev/null)
MACROS=$(addprefix -D,SVN_VERS=\"$(SVN_VERS)\")
INCLUDE_PATH=src
LIBRARY_PATH=

GCC=gcc

################################################ define the projects ############################################

#collect all the projects
projects=$(addprefix projects/, \
	bubbles/tm_disconnected \
	eight_BK/smeared_BK_all_in_one \
	g/g \
	generate_confs/generate_stag \
	nucleons/nucleons_smeared \
	reno_const/RIMOM \
	semileptonic/new_semileptonic \
	semileptonic/semileptonic_smeared \
	semileptonic/semileptonic_smeared_point_source \
	static_potential/compute_potential)

tools=$(addprefix tools/, \
	ave \
	clusterize2 \
	endianess_check/endianess_check \
	print_gamma/gamma_test \
	print_gamma/gamma_comm \
	unitarity_check/unitarity_check \
	spincolor_compare/spincolor_compare \
	meson_2pts/meson_2pts \
	meson_2pts/meson_2pts_point_source \
	conf_convert/ildg_to_eo \
	conf_convert/generate_random)

############################################# define the library pieces #########################################

#add single files into library pieces
base=$(addprefix base/, \
	close communicate debug \
	global_variables \
	init \
	random \
	routines \
	vectors)
dirac_operators=$(addprefix dirac_operators/dirac_operator_, \
	stD/dirac_operator_stD \
	tmDeoimpr/dirac_operator_tmDeoimpr \
	tmQ/dirac_operator_tmQ \
	tmQ_left/dirac_operator_tmQ_left \
	tmQ2/dirac_operator_tmQ2 \
	tmQ/dirac_operator_tmQ_128 \
	tmQ/reconstruct_tm_doublet \
	tmQ2/dirac_operator_tmQ2_128)
geometry=$(addprefix geometry/, \
	geometry_eo \
	geometry_lx \
	geometry_mix)
hmc=$(addprefix hmc/, \
	backfield \
	gauge/wilson_force \
	momenta/momenta_generation \
	momenta/momenta_action \
	rootst_eoimpr/rootst_eoimpr_action \
	rootst_eoimpr/rootst_eoimpr_eigenvalues \
	rootst_eoimpr/rootst_eoimpr_force \
	rootst_eoimpr/rootst_eoimpr_omelyan_integrator \
	rootst_eoimpr/rootst_eoimpr_pseudofermions_generation \
	rootst_eoimpr/rootst_eoimpr_rhmc_step)
inverters=$(addprefix inverters/twisted_mass/, \
	cg_invert_tmDeoimpr \
	cg_invert_tmQ2 \
	cg_128_invert_tmQ2 \
	cgm_invert_tmQ2 tm_frontends) \
	$(addprefix inverters/staggered/, \
	cgm_invert_stD2ee_m2 \
	cg_invert_stD2ee_m2 \
	cg_invert_stD2_m2 \
	cg_invert_evn_stD)
IO=$(addprefix IO/, \
	checksum \
	endianess \
	ILDG_File \
	input \
	reader \
	writer)
new_types=$(addprefix new_types/, \
	complex \
	dirac \
	float128 \
	rat_exp \
	spin \
	su3)
operations=$(addprefix operations/, \
	contract \
	fft \
	fourier_transform \
	gauge_fixing \
	gaugeconf \
	remap_vector \
	smear \
	su3_paths \
	vector_gather)

#collect all pieces of the library and define targets for library pieces
nissa_library_pieces=$(addprefix src/, $(base) $(dirac_operators) $(geometry) $(hmc) $(inverters) $(IO) $(new_types) $(operations) linalgs/linalgs)

################################################## global targets ###############################################

all: Makefile src/libnissa.a $(projects) $(tools)

Makefile: default.Makefile configure 
	./configure; \
	echo "WARNING: Makefile update, rerun it!!!!"; \
	exit 1

#define targets for library pieces
nissa_library_objects: $(addsuffix .o,$(nissa_library_pieces)) Makefile

#include table of dependencies
-include $(addsuffix .d,$(nissa_library_pieces) $(projects) $(tools))

################################## rules to produce objects, library and projects ################################

$(projects) $(tools): %: %.o src/libnissa.a Makefile
	$(CC)                           \
	$(addprefix -I,$(INCLUDE_PATH)) \
	$(addprefix -L,$(LIBRARY_PATH)) \
	$(CFLAGS)                       \
	$(MACROS)                       \
	-o $@                           \
	$(addsuffix .o, $@)             \
	src/libnissa.a

src/libnissa.a: $(addsuffix .o,$(nissa_library_pieces)) Makefile
	ar cru src/libnissa.a $(addsuffix .o,$(nissa_library_pieces))

$(addsuffix .o,$(nissa_library_pieces) $(projects) $(tools)): %.o: %.cpp %.d Makefile
	$(CC)				\
	$(addprefix -I,$(INCLUDE_PATH)) \
	$<				\
	$(CFLAGS) 			\
	$(MACROS) 			\
	-c -o $@

$(addsuffix .d,$(nissa_library_pieces) $(projects) $(tools)): %.d: %.cpp Makefile
	$(GCC) $< -MM -MT $(@:%.d=%.o) $(addprefix -I,$(INCLUDE_PATH)) -o $@


clean:
	 rm -rf $(addsuffix .d,$(nissa_library_pieces)) $(addsuffix .o,$(nissa_library_pieces)) src/libnissa.a $(projects) $(tools)

.PHONY: clean all nissa_library_objects
