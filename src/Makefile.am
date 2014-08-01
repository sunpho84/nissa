include_HEADERS=nissa.hpp

#add nissa library
lib_LIBRARIES=$(top_builddir)/lib/libnissa.a

########################################### base sources ##################################

__top_builddir__lib_libnissa_a_SOURCES= \
        base/bench.cpp \
        base/close.cpp \
        base/debug.cpp \
        base/global_variables.cpp \
        base/init.cpp \
        base/random.cpp \
        base/vectors.cpp
include_HEADERS+= \
        base/bench.hpp \
        base/close.hpp \
        base/debug.hpp \
        base/global_variables.hpp \
        base/init.hpp \
        base/macros.hpp \
        base/random.hpp \
        base/thread_macros.hpp \
        base/vectors.hpp

#add BGQ if needed
if BGQ
include_HEADERS+=base/bgq_instructions.hpp
endif

#svn versioning: specify to build sources and include script
BUILT_SOURCES=update_svnversion svnversion.hpp
update_svnversion:
	$(srcdir)/base/update_svnversion.hpp.sh $(top_srcdir)
.PHONY: update_svnversion
EXTRA_DIST=svnversion.hpp base/update_svnversion.hpp.sh

#include global_variables.cpp in the installation
#needed because included in global_variables.hpp
global_variabledir=$(prefix)/include/base
global_variable_DATA=base/global_variables.cpp

#clean svnversion.hpp
CLEANFILES=base/svnversion.hpp

######################################### bgq specific #####################################
if BGQ
__top_builddir__lib_libnissa_a_SOURCES+= \
        bgq/Wilson_hopping_matrix_lx_bgq.cpp \
        bgq/staggered_hopping_matrix_eo_or_oe_bgq.cpp

include_HEADERS+= \
	bgq/Wilson_hopping_matrix_lx_bgq.hpp \
        bgq/staggered_hopping_matrix_eo_or_oe_bgq.hpp \
        bgq/bgq_macros.hpp \
        bgq/intrinsic/declare.hpp \
        bgq/intrinsic/load.hpp \
        bgq/intrinsic/mergesplit.hpp \
        bgq/intrinsic/oper.hpp \
        bgq/intrinsic/store.hpp \
        bgq/intrinsic.hpp

#bgq barrier
if BGQ_EMU
else
__top_builddir__lib_libnissa_a_SOURCES+=bgq/bgq_barrier.c
include_HEADERS+=bgq/bgq_barrier.hpp
endif

#compile or not spi stuff
if SPI
__top_builddir__lib_libnissa_a_SOURCES+=bgq/spi.cpp
include_HEADERS+=bgq/spi.hpp
endif

endif

######################################### communicate ######################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	communicate/all_to_all.cpp \
	communicate/borders.cpp \
	communicate/edges.cpp

include_HEADERS+= \
	communicate/all_to_all.hpp \
	communicate/borders.hpp \
	communicate/communicate.hpp \
	communicate/edges.hpp

######################################## dirac operators#####################################

__top_builddir__lib_libnissa_a_SOURCES+= \
        dirac_operators/stD/dirac_operator_stD.cpp \
        dirac_operators/tmDeoimpr/dirac_operator_tmDeoimpr.cpp \
        dirac_operators/tmQ/dirac_operator_tmQ.cpp \
        dirac_operators/tmQ_left/dirac_operator_tmQ_left.cpp \
        dirac_operators/tmQ2/dirac_operator_tmQ2.cpp \
        dirac_operators/tmQ/reconstruct_tm_doublet.cpp \
        dirac_operators/tmclovQ/reconstruct_tmclov_doublet.cpp \
        dirac_operators/Wstat/dirac_operator_Wstat.cpp \
        dirac_operators/WclovQ/dirac_operator_WclovQ.cpp \
        dirac_operators/WclovQ2/dirac_operator_WclovQ2.cpp \
        dirac_operators/tmclovQ/dirac_operator_tmclovQ.cpp \
        dirac_operators/tmclovQ2/dirac_operator_tmclovQ2.cpp \
        dirac_operators/tmQ/dirac_operator_tmQ_128.cpp \
        dirac_operators/tmDeoimpr/dirac_operator_tmDeoimpr_128.cpp \
        dirac_operators/tmQ2/dirac_operator_tmQ2_128.cpp

include_HEADERS+= \
	dirac_operators/stD/dirac_operator_stD.hpp\
	dirac_operators/tmDeoimpr/dirac_operator_tmDeoimpr.hpp \
        dirac_operators/tmDeoimpr/dirac_operator_tmDeoimpr_128.hpp \
	dirac_operators/tmQ/dirac_operator_tmQ.hpp \
        dirac_operators/tmQ/reconstruct_tm_doublet.hpp \
        dirac_operators/tmQ/dirac_operator_tmQ_128.hpp \
	dirac_operators/tmQ_left/dirac_operator_tmQ_left.hpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2.hpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2_128.hpp \
	dirac_operators/tmclovQ/dirac_operator_tmclovQ.hpp \
	dirac_operators/tmclovQ/reconstruct_tmclov_doublet.hpp \
	dirac_operators/Wstat/dirac_operator_Wstat.hpp \
	dirac_operators/WclovQ/dirac_operator_WclovQ.hpp \
	dirac_operators/WclovQ2/dirac_operator_WclovQ2.hpp \
	dirac_operators/tmclovQ2/dirac_operator_tmclovQ2.hpp

#compile or not bgq source
if BGQ
__top_builddir__lib_libnissa_a_SOURCES+= \
        dirac_operators/tmQ/dirac_operator_tmQ_bgq.cpp \
        dirac_operators/tmQ2/dirac_operator_tmQ2_bgq.cpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2_128_bgq.cpp \
	dirac_operators/stD/dirac_operator_stD_bgq.cpp
include_HEADERS+= \
        dirac_operators/tmQ/dirac_operator_tmQ_bgq.hpp \
        dirac_operators/tmQ2/dirac_operator_tmQ2_bgq.hpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2_128_bgq.hpp \
	dirac_operators/stD/dirac_operator_stD_bgq.hpp
endif

EXTRA_DIST+= \
        dirac_operators/WclovQ/dirac_operator_WclovQ_portable.cpp \
        dirac_operators/stD/dirac_operator_stD_portable.cpp \
        dirac_operators/stD/dirac_operator_stD_32_portable.cpp \
        dirac_operators/tmDeoimpr/dirac_operator_tmDeoimpr_portable.cpp \
	dirac_operators/tmQ/dirac_operator_tmQ_portable.cpp \
        dirac_operators/tmQ_left/dirac_operator_tmQ_left_portable.cpp

################################### geometry ###################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	geometry/geometry_eo.cpp \
	geometry/geometry_lx.cpp \
	geometry/geometry_mix.cpp \
	geometry/geometry_Wsklx.cpp

include_HEADERS+= \
	geometry/geometry_eo.hpp \
	geometry/geometry_lx.hpp \
	geometry/geometry_mix.hpp \
	geometry/geometry_Wsklx.hpp

#compile or not virtual parallelization
if USE_VNODES
__top_builddir__lib_libnissa_a_SOURCES+=geometry/geometry_vir.cpp
include_HEADERS+=geometry/geometry_vir.hpp
endif

###################################### gpu #########################################

if GPU

include $(top_srcdir)/config/cuda.mk

LINK=$(LIBTOOL) --mode=link $(CC) -o $@ $(CUDA_LDFLAGS) $(CUDA_LIBS)

__top_builddir__lib_libnissa_a_SOURCES+= \
	gpu/fields.cu \
	gpu/global_variables.cu \
	gpu/init.cu \
	gpu/new_types.cu \
	gpu/wrappers.cu

include_HEADERS+= \
	gpu/fields.hpp \
	gpu/global_variables.hpp \
	gpu/init.hpp \
	gpu/macros.hpp \
	gpu/new_types.hpp \
        gpu/wrappers.hpp

endif

###################################### hmc #########################################

__top_builddir__lib_libnissa_a_SOURCES+= \
        hmc/backfield.cpp \
        hmc/gauge/gluonic_force.cpp \
        hmc/gauge/topological_action.cpp \
        hmc/gauge/topological_force.cpp \
        hmc/gauge/tree_level_Symanzik_force.cpp \
        hmc/gauge/tree_level_Symanzik_action.cpp \
        hmc/gauge/Wilson_force.cpp \
        hmc/momenta/momenta_generation.cpp \
        hmc/momenta/momenta_action.cpp \
        hmc/rootst_eoimpr/rootst_eoimpr_action.cpp \
        hmc/rootst_eoimpr/rootst_eoimpr_eigenvalues.cpp \
        hmc/rootst_eoimpr/rootst_eoimpr_quark_force.cpp \
        hmc/rootst_eoimpr/rootst_eoimpr_omelyan_integrator.cpp \
        hmc/rootst_eoimpr/rootst_eoimpr_pseudofermions_generation.cpp \
        hmc/rootst_eoimpr/rootst_eoimpr_rhmc_step.cpp

include_HEADERS+= \
        hmc/backfield.hpp \
        hmc/gauge/gluonic_force.hpp \
        hmc/gauge/topological_action.hpp \
        hmc/gauge/topological_force.hpp \
        hmc/gauge/tree_level_Symanzik_force.hpp \
        hmc/gauge/tree_level_Symanzik_action.hpp \
        hmc/gauge/Wilson_force.hpp \
        hmc/momenta/momenta_generation.hpp \
        hmc/momenta/momenta_action.hpp \
        hmc/rootst_eoimpr/rootst_eoimpr_action.hpp \
        hmc/rootst_eoimpr/rootst_eoimpr_eigenvalues.hpp \
        hmc/rootst_eoimpr/rootst_eoimpr_quark_force.hpp \
        hmc/rootst_eoimpr/rootst_eoimpr_omelyan_integrator.hpp \
        hmc/rootst_eoimpr/rootst_eoimpr_pseudofermions_generation.hpp \
        hmc/rootst_eoimpr/rootst_eoimpr_rhmc_step.hpp

EXTRA_DIST+=hmc/rootst_eoimpr/rat_expansion_database.cpp

###################################### inverters #################################

__top_builddir__lib_libnissa_a_SOURCES+= \
        inverters/twisted_mass/cg_invert_tmDeoimpr.cpp \
        inverters/twisted_mass/cg_invert_tmQ2.cpp \
        inverters/twisted_mass/cg_64_invert_tmQ2.cpp \
        inverters/twisted_mass/cgm_invert_tmQ2.cpp \
        inverters/twisted_mass/cg_128_invert_tmQ2.cpp \
        inverters/twisted_mass/cg_128_invert_tmDeoimpr.cpp \
        inverters/twisted_mass/tm_frontends.cpp \
        inverters/staggered/cgm_invert_stD2ee_m2_portable.cpp \
        inverters/staggered/cgm_32_invert_stD2ee_m2_portable.cpp \
        inverters/staggered/cg_invert_stD2ee_m2_portable.cpp \
        inverters/staggered/cgm_invert_stD2ee_m2.cpp \
        inverters/staggered/cg_invert_stD2ee_m2.cpp \
        inverters/staggered/cg_invert_stD.cpp \
        inverters/staggered/cg_invert_evn_stD.cpp \
        inverters/Wclov/cg_invert_WclovQ.cpp \
        inverters/Wclov/cg_invert_WclovQ2.cpp \
        inverters/tmclov/cg_invert_tmclovQ.cpp \
        inverters/tmclov/cg_invert_tmclovQ2.cpp \
        inverters/tmclov/cgm_invert_tmclovQ2.cpp \
        inverters/tmclov/tmclov_frontends.cpp

include_HEADERS+= \
        inverters/twisted_mass/cg_invert_tmDeoimpr.hpp \
        inverters/twisted_mass/cg_invert_tmQ2.hpp \
        inverters/twisted_mass/cg_64_invert_tmQ2.hpp \
        inverters/twisted_mass/cgm_invert_tmQ2.hpp \
        inverters/twisted_mass/cg_128_invert_tmQ2.hpp \
        inverters/twisted_mass/cg_128_invert_tmDeoimpr.hpp \
        inverters/twisted_mass/tm_frontends.hpp \
        inverters/staggered/cgm_invert_stD2ee_m2.hpp \
        inverters/staggered/cg_invert_stD2ee_m2.hpp \
        inverters/staggered/cg_invert_stD.hpp \
        inverters/staggered/cg_invert_evn_stD.hpp \
	inverters/staggered/cgm_invert_stD2ee_m2_portable.hpp \
	inverters/staggered/cg_invert_stD2ee_m2_portable.hpp \
        inverters/Wclov/cg_invert_WclovQ.hpp \
        inverters/Wclov/cg_invert_WclovQ2.hpp \
        inverters/tmclov/cg_invert_tmclovQ.hpp \
        inverters/tmclov/cg_invert_tmclovQ2.hpp \
        inverters/tmclov/cgm_invert_tmclovQ2.hpp \
        inverters/tmclov/tmclov_frontends.hpp

#compile or not bgq source
if BGQ
__top_builddir__lib_libnissa_a_SOURCES+= \
        inverters/staggered/cgm_invert_stD2ee_m2_bgq.cpp \
        inverters/staggered/cg_invert_stD2ee_m2_bgq.cpp \
        inverters/twisted_mass/cgm_invert_tmQ2_bgq.cpp \
        inverters/twisted_mass/cg_64_invert_tmQ2_bgq.cpp \
        inverters/twisted_mass/cg_128_invert_tmQ2_bgq.cpp \
	inverters/twisted_mass/cg_invert_tmQ2_bgq.cpp

include_HEADERS+= \
	inverters/staggered/cgm_invert_stD2ee_m2_bgq.hpp \
	inverters/staggered/cg_invert_stD2ee_m2_bgq.hpp \
        inverters/twisted_mass/cgm_invert_tmQ2_bgq.hpp \
        inverters/twisted_mass/cg_64_invert_tmQ2_bgq.hpp \
        inverters/twisted_mass/cg_128_invert_tmQ2_bgq.hpp \
	inverters/twisted_mass/cg_invert_tmQ2_bgq.hpp
endif

EXTRA_DIST+= \
        inverters/twisted_mass/cg_invert_tmQ2_common.cpp \
	inverters/twisted_mass/cgm_invert_tmQ2_bgq.hpp \
        inverters/templates/cg_128_invert_template_threaded.cpp \
        inverters/templates/cg_invert_template_threaded.cpp \
        inverters/templates/cgm_invert_template_threaded.cpp

############################################### io #################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	io/ILDG_File.cpp \
	io/checksum.cpp \
	io/endianness.cpp \
	io/input.cpp \
	io/reader.cpp \
	io/writer.cpp

include_HEADERS+= \
	io/ILDG_File.hpp \
	io/checksum.hpp \
	io/endianness.hpp \
	io/input.hpp \
	io/reader.hpp \
	io/writer.hpp

########################################### linalgs ##################################

__top_builddir__lib_libnissa_a_SOURCES+=linalgs/linalgs.cpp
include_HEADERS+=linalgs/linalgs.hpp

########################################## new types #################################

__top_builddir__lib_libnissa_a_SOURCES+= \
        new_types/complex.cpp \
        new_types/dirac.cpp \
        new_types/float_128.cpp \
        new_types/float_256.cpp \
        new_types/rat_approx.cpp \
        new_types/read_new_types.cpp \
        new_types/spin.cpp \
        new_types/su3.cpp

include_HEADERS+= \
	new_types/complex.hpp \
        new_types/dirac.hpp \
        new_types/float_128.hpp \
        new_types/float_256.hpp \
	new_types/new_types_definitions.hpp \
        new_types/rat_approx.hpp \
        new_types/read_new_types.hpp \
        new_types/spin.hpp \
        new_types/su3.hpp

####################################### operations ##################################

__top_builddir__lib_libnissa_a_SOURCES+= \
        operations/contract/mesons_2pts.cpp \
        operations/contract/mesons_eight.cpp \
        operations/contract/optimized_mesons_2pts.cpp \
        operations/contract/site_contract.cpp \
        operations/contract/stag.cpp \
        operations/covariant_derivative.cpp \
        operations/fft.cpp \
        operations/fourier_transform.cpp \
        operations/gauge_fixing.cpp \
        operations/gaugeconf.cpp \
        operations/remap_vector.cpp \
        operations/remez/remez_algorithm.cpp \
        operations/shift.cpp \
        operations/source.cpp \
        operations/smearing/APE.cpp \
        operations/smearing/HYP.cpp \
        operations/smearing/stout_common.cpp \
        operations/smearing/stout_eo.cpp \
        operations/smearing/stout_lx.cpp \
        operations/smearing/gaussian.cpp \
        operations/su3_paths/all_rectangles.cpp \
        operations/su3_paths/arbitrary.cpp \
	operations/su3_paths/gauge_sweeper.cpp \
	operations/su3_paths/hypercubic_transfer.cpp \
        operations/su3_paths/plaquette.cpp \
        operations/su3_paths/pline.cpp \
        operations/su3_paths/rectangles.cpp \
        operations/su3_paths/rectangular_staples.cpp \
        operations/su3_paths/shortest_hypercubic_paths.cpp \
        operations/su3_paths/squared_staples.cpp \
        operations/su3_paths/topological_charge.cpp \
        operations/vector_gather.cpp

include_HEADERS+= \
        operations/covariant_derivative.hpp \
        operations/fft.hpp \
        operations/fourier_transform.hpp \
        operations/gauge_fixing.hpp \
        operations/gaugeconf.hpp \
        operations/remap_vector.hpp \
        operations/shift.hpp \
        operations/source.hpp \
        operations/vector_gather.hpp \
	operations/remez/remez_algorithm.hpp \
        operations/contract/mesons_2pts.hpp \
        operations/contract/mesons_eight.hpp \
        operations/contract/optimized_mesons_2pts.hpp \
        operations/contract/site_contract.hpp \
        operations/contract/stag.hpp \
        operations/smearing/APE.hpp \
        operations/smearing/HYP.hpp \
        operations/smearing/stout.hpp \
        operations/smearing/gaussian.hpp \
        operations/su3_paths/all_rectangles.hpp \
        operations/su3_paths/arbitrary.hpp \
	operations/su3_paths/gauge_sweeper.hpp \
	operations/su3_paths/hypercubic_transfer.hpp \
        operations/su3_paths/plaquette.hpp \
        operations/su3_paths/pline.hpp \
        operations/su3_paths/rectangles.hpp \
        operations/su3_paths/rectangular_staples.hpp \
        operations/su3_paths/shortest_hypercubic_paths.hpp \
        operations/su3_paths/squared_staples.hpp \
        operations/su3_paths/topological_charge.hpp

########################################## routines ###################################

__top_builddir__lib_libnissa_a_SOURCES+= \
        routines/ios.cpp \
        routines/math_routines.cpp \
        routines/mpi_routines.cpp

include_HEADERS+= \
	routines/ios.hpp \
        routines/math_routines.hpp \
        routines/mpi_routines.hpp

#compile or not thread support
if USE_THREADS
__top_builddir__lib_libnissa_a_SOURCES+=routines/thread.cpp
include_HEADERS+=routines/thread.hpp
endif
