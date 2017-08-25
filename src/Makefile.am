include_HEADERS=nissa.hpp

#add nissa library
lib_LIBRARIES=$(top_builddir)/lib/libnissa.a

base/git_info.hpp: $(top_srcdir)/.git/HEAD $(top_srcdir)/.git/index $(top_srcdir)/config/get_git_info.sh Makefile
	@echo "Reading git version"
	@$(top_srcdir)/config/get_git_info.sh $(top_srcdir) > $@
$(top_srcdir)/.git/HEAD:
$(top_srcdir)/.git/index:


BUILT_SOURCES=base/git_info.hpp

########################################### base sources ##################################

__top_builddir__lib_libnissa_a_SOURCES= \
	base/bench.cpp \
	base/close.cpp \
	base/debug.cpp \
	base/init.cpp \
	base/random.cpp \
	base/vectors.cpp
include_HEADERS+= \
	base/bench.hpp \
	base/close.hpp \
	base/debug.hpp \
	base/git_info.hpp \
	base/init.hpp \
	base/random.hpp \
	base/thread_macros.hpp \
	base/vectors.hpp

if USE_TMLQCD
__top_builddir__lib_libnissa_a_SOURCES+= \
	base/tmLQCD_bridge.cpp
include_HEADERS+= \
	base/tmLQCD_bridge.hpp
endif

if USE_DDALPHAAMG
__top_builddir__lib_libnissa_a_SOURCES+= \
	base/DDalphaAMG_bridge.cpp
include_HEADERS+= \
	base/DDalphaAMG_bridge.hpp
endif

#add BGQ if needed
if BGQ
include_HEADERS+=base/bgq_instructions.hpp
endif

######################################### bgq specific #####################################
if BGQ
__top_builddir__lib_libnissa_a_SOURCES+= \
	bgq/Wilson_hopping_matrix_lx_bgq.cpp \
	bgq/Wilson_hopping_matrix_eo_or_oe_bgq.cpp \
	bgq/staggered_hopping_matrix_eo_or_oe_bgq.cpp

include_HEADERS+= \
	bgq/Wilson_hopping_matrix_lx_bgq.hpp \
	bgq/Wilson_hopping_matrix_eo_or_oe_bgq_32_64.hpp \
	bgq/Wilson_hopping_matrix_eo_or_oe_bgq.hpp \
	bgq/staggered_hopping_matrix_eo_or_oe_bgq_32_64.hpp \
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
	communicate/communicate.cpp \
	communicate/edges.cpp

include_HEADERS+= \
	communicate/all_to_all.hpp \
	communicate/borders.hpp \
	communicate/communicate.hpp \
	communicate/edges.hpp

######################################## dirac operators#####################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	dirac_operators/momenta/MFACC.cpp \
	dirac_operators/stD/dirac_operator_stD.cpp \
	dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec.cpp \
	dirac_operators/tmQ/dirac_operator_tmQ.cpp \
	dirac_operators/tmQ_left/dirac_operator_tmQ_left.cpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2.cpp \
	dirac_operators/tmQ/reconstruct_tm_doublet.cpp \
	dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec.cpp \
	dirac_operators/tmclovQ/reconstruct_tmclov_doublet.cpp \
	dirac_operators/Wstat/dirac_operator_Wstat.cpp \
	dirac_operators/WclovQ/dirac_operator_WclovQ.cpp \
	dirac_operators/WclovQ2/dirac_operator_WclovQ2.cpp \
	dirac_operators/tmclovQ/dirac_operator_tmclovQ.cpp \
	dirac_operators/tmclovQ/dirac_operator_tmclovQ_128.cpp \
	dirac_operators/tmclovQ2/dirac_operator_tmclovQ2.cpp \
	dirac_operators/tmclovQ2/dirac_operator_tmclovQ2_128.cpp \
	dirac_operators/tmQ/dirac_operator_tmQ_128.cpp \
	dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_128.cpp \
	dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec_128.cpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2_128.cpp

include_HEADERS+= \
	dirac_operators/momenta/MFACC.hpp \
	dirac_operators/stD/dirac_operator_stD.hpp\
	dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec.hpp \
	dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_128.hpp \
	dirac_operators/tmQ/dirac_operator_tmQ.hpp \
	dirac_operators/tmQ/reconstruct_tm_doublet.hpp \
	dirac_operators/tmQ/dirac_operator_tmQ_128.hpp \
	dirac_operators/tmQ_left/dirac_operator_tmQ_left.hpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2.hpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2_128.hpp \
	dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec.hpp \
	dirac_operators/tmclovQ/dirac_operator_tmclovQ.hpp \
	dirac_operators/tmclovQ/dirac_operator_tmclovQ_128.hpp \
	dirac_operators/tmclovQ/reconstruct_tmclov_doublet.hpp \
	dirac_operators/Wstat/dirac_operator_Wstat.hpp \
	dirac_operators/WclovQ/dirac_operator_WclovQ.hpp \
	dirac_operators/WclovQ2/dirac_operator_WclovQ2.hpp \
	dirac_operators/tmclovQ2/dirac_operator_tmclovQ2.hpp \
	dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec_128.hpp \
	dirac_operators/tmclovQ2/dirac_operator_tmclovQ2_128.hpp

#compile or not bgq source
if BGQ
__top_builddir__lib_libnissa_a_SOURCES+= \
	dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_bgq.cpp \
	dirac_operators/tmQ/dirac_operator_tmQ_bgq.cpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2_bgq.cpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2_128_bgq.cpp \
	dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec_bgq.cpp \
	dirac_operators/tmclovQ/dirac_operator_tmclovQ_bgq.cpp \
	dirac_operators/tmclovQ2/dirac_operator_tmclovQ2_bgq.cpp \
	dirac_operators/tmclovQ2/dirac_operator_tmclovQ2_128_bgq.cpp \
	dirac_operators/stD/dirac_operator_stD_bgq.cpp
include_HEADERS+= \
	dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_bgq.hpp \
	dirac_operators/tmQ/dirac_operator_tmQ_bgq.hpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2_bgq.hpp \
	dirac_operators/tmQ2/dirac_operator_tmQ2_128_bgq.hpp \
	dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec_bgq.hpp \
	dirac_operators/tmclovQ/dirac_operator_tmclovQ_bgq.hpp \
	dirac_operators/tmclovQ2/dirac_operator_tmclovQ2_bgq.hpp \
	dirac_operators/tmclovQ2/dirac_operator_tmclovQ2_128_bgq.hpp \
	dirac_operators/stD/dirac_operator_stD_bgq.hpp
endif

EXTRA_DIST= \
	base/git_info.hpp \
	dirac_operators/WclovQ/dirac_operator_WclovQ_portable.cpp \
	dirac_operators/stD/dirac_operator_stD_portable.cpp \
	dirac_operators/stD/dirac_operator_stD_32_portable.cpp \
	dirac_operators/stD/dirac_operator_stD_bgq_template.cpp \
	dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec_portable.cpp \
	dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_portable.cpp \
	dirac_operators/tmQ/dirac_operator_tmQ_portable.cpp \
	dirac_operators/tmQ_left/dirac_operator_tmQ_left_portable.cpp

################################## free theory##################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	free_theory/cg_eoprec_twisted_free_operator.cpp \
	free_theory/free_theory_types_routines.cpp \
	free_theory/tlSym_gauge_propagator.cpp \
	free_theory/twisted_free_Dirac_eoprec_operator.cpp \
	free_theory/twisted_propagator.cpp

include_HEADERS+= \
	free_theory/cg_eoprec_twisted_free_operator.hpp \
	free_theory/free_theory_types.hpp \
	free_theory/free_theory_types_routines.hpp \
	free_theory/tlSym_gauge_propagator.hpp \
	free_theory/twisted_free_Dirac_eoprec_operator.hpp \
	free_theory/twisted_propagator.hpp

################################### geometry ###################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	geometry/geometry_eo.cpp \
	geometry/geometry_Leb.cpp \
	geometry/geometry_lx.cpp \
	geometry/geometry_mix.cpp

include_HEADERS+= \
	geometry/geometry_eo.hpp \
	geometry/geometry_Leb.hpp \
	geometry/geometry_lx.hpp \
	geometry/geometry_mix.hpp

#compile or not virtual parallelization
if USE_VNODES
__top_builddir__lib_libnissa_a_SOURCES+=geometry/geometry_vir.cpp
include_HEADERS+=geometry/geometry_vir.hpp
endif

###################################### hmc #########################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	hmc/backfield.cpp \
	hmc/theory_pars.cpp \
	hmc/fermions/rootst_eoimpr_quark_force.cpp \
	hmc/fermions/pseudofermions_generation.cpp \
	hmc/gauge/MFACC_fields.cpp \
	hmc/gauge/gluonic_force.cpp \
	hmc/gauge/pure_gauge_hmc_step.cpp \
	hmc/gauge/pure_gauge_implicit_integrator.cpp \
	hmc/gauge/pure_gauge_Omelyan_integrator.cpp \
	hmc/gauge/Symanzik_force.cpp \
	hmc/gauge/Symanzik_action.cpp \
	hmc/gauge/topological_action.cpp \
	hmc/gauge/topological_force.cpp \
	hmc/gauge/Wilson_action.cpp \
	hmc/gauge/Wilson_force.cpp \
	hmc/momenta/momenta_action.cpp \
	hmc/momenta/momenta_evolve.cpp \
	hmc/momenta/momenta_generation.cpp \
	hmc/multipseudo/multipseudo_rhmc_step.cpp \
	hmc/multipseudo/Omelyan_integrator.cpp \
	hmc/multipseudo/quark_force.cpp \
	hmc/multipseudo/set_expansions.cpp \
	hmc/multipseudo/theory_action.cpp

include_HEADERS+= \
	hmc/backfield.hpp \
	hmc/hmc.hpp \
	hmc/quark_pars.hpp \
	hmc/theory_pars.hpp \
	hmc/fermions/rootst_eoimpr_quark_force.hpp \
	hmc/fermions/pseudofermions_generation.hpp \
	hmc/gauge/MFACC_fields.hpp \
	hmc/gauge/gluonic_action.hpp \
	hmc/gauge/gluonic_force.hpp \
	hmc/gauge/pure_gauge_hmc_step.hpp \
	hmc/gauge/pure_gauge_implicit_integrator.hpp \
	hmc/gauge/pure_gauge_Omelyan_integrator.hpp \
	hmc/gauge/Symanzik_force.hpp \
	hmc/gauge/Symanzik_action.hpp \
	hmc/gauge/topological_action.hpp \
	hmc/gauge/topological_force.hpp \
	hmc/gauge/Wilson_action.hpp \
	hmc/gauge/Wilson_force.hpp \
	hmc/momenta/momenta_action.hpp \
	hmc/momenta/momenta_evolve.hpp \
	hmc/momenta/momenta_generation.hpp \
	hmc/multipseudo/multipseudo_rhmc_step.hpp \
	hmc/multipseudo/Omelyan_integrator.hpp \
	hmc/multipseudo/quark_force.hpp \
	hmc/multipseudo/set_expansions.hpp \
	hmc/multipseudo/theory_action.hpp

###################################### inverters #################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	inverters/momenta/cg_invert_MFACC.cpp \
	inverters/momenta/cgm_invert_MFACC.cpp \
	inverters/staggered/cgm_invert_stD2ee_m2_portable.cpp \
	inverters/staggered/cgm_32_invert_stD2ee_m2.cpp \
	inverters/staggered/cgm_32_invert_stD2ee_m2_portable.cpp \
	inverters/staggered/cg_invert_stD2ee_m2_portable.cpp \
	inverters/staggered/cg_invert_stD2Leb_ee_m2_portable.cpp \
	inverters/staggered/cgm_invert_stD2ee_m2.cpp \
	inverters/staggered/cg_invert_stD2ee_m2.cpp \
	inverters/staggered/cg_invert_stD.cpp \
	inverters/staggered/cg_invert_evn_stD.cpp \
	inverters/twisted_mass/cg_invert_tmD_eoprec.cpp \
	inverters/twisted_mass/cg_invert_tmQ2.cpp \
	inverters/twisted_mass/cg_64_invert_tmQ2.cpp \
	inverters/twisted_mass/cg_64_invert_tmD_eoprec.cpp \
	inverters/twisted_mass/cgm_invert_tmQ2.cpp \
	inverters/twisted_mass/cg_128_invert_tmQ2.cpp \
	inverters/twisted_mass/cg_128_invert_tmD_eoprec.cpp \
	inverters/twisted_mass/tm_frontends.cpp \
	inverters/twisted_clover/cg_invert_tmclovQ2.cpp \
	inverters/twisted_clover/cg_64_invert_tmclovQ2_portable.cpp \
	inverters/twisted_clover/cg_64_invert_tmclovD_eoprec.cpp \
	inverters/twisted_clover/cg_128_invert_tmclovQ2.cpp \
	inverters/twisted_clover/cg_128_invert_tmclovD_eoprec.cpp \
	inverters/twisted_clover/cgm_invert_tmclovDkern_eoprec_square_portable.cpp \
	inverters/twisted_clover/cgm_invert_tmclovQ2.cpp \
	inverters/twisted_clover/tmclov_frontends.cpp \
	inverters/Wclov/cg_invert_WclovQ.cpp \
	inverters/Wclov/cg_invert_WclovQ2.cpp \
	inverters/twisted_clover/cg_invert_tmclovD_eoprec.cpp \
	inverters/twisted_clover/cg_invert_tmclovQ.cpp

include_HEADERS+= \
	inverters/momenta/cg_invert_MFACC.hpp \
	inverters/momenta/cgm_invert_MFACC.hpp \
	inverters/twisted_mass/cg_invert_tmD_eoprec.hpp \
	inverters/twisted_mass/cg_invert_tmQ2.hpp \
	inverters/twisted_mass/cg_64_invert_tmQ2.hpp \
	inverters/twisted_mass/cg_64_invert_tmD_eoprec.hpp \
	inverters/twisted_mass/cgm_invert_tmQ2.hpp \
	inverters/twisted_mass/cg_128_invert_tmQ2.hpp \
	inverters/twisted_mass/cg_128_invert_tmD_eoprec.hpp \
	inverters/twisted_mass/tm_frontends.hpp \
	inverters/staggered/cgm_invert_stD2ee_m2.hpp \
	inverters/staggered/cg_invert_stD2ee_m2.hpp \
	inverters/staggered/cg_invert_stD.hpp \
	inverters/staggered/cg_invert_evn_stD.hpp \
	inverters/staggered/cgm_32_invert_stD2ee_m2.hpp \
	inverters/staggered/cgm_invert_stD2ee_m2_portable.hpp \
	inverters/staggered/cgm_32_invert_stD2ee_m2_portable.hpp \
	inverters/staggered/cg_invert_stD2ee_m2_portable.hpp \
	inverters/staggered/cg_invert_stD2Leb_ee_m2_portable.hpp \
	inverters/Wclov/cg_invert_WclovQ.hpp \
	inverters/Wclov/cg_invert_WclovQ2.hpp \
	inverters/twisted_clover/cg_invert_tmclovQ.hpp \
	inverters/twisted_clover/cg_invert_tmclovQ2.hpp \
	inverters/twisted_clover/cg_64_invert_tmclovQ2_portable.hpp \
	inverters/twisted_clover/cg_64_invert_tmclovQ2.hpp \
	inverters/twisted_clover/cg_64_invert_tmclovD_eoprec.hpp \
	inverters/twisted_clover/cg_128_invert_tmclovQ2.hpp \
	inverters/twisted_clover/cg_128_invert_tmclovD_eoprec.hpp \
	inverters/twisted_clover/cgm_invert_tmclovQ2.hpp \
	inverters/twisted_clover/cg_invert_tmclovD_eoprec.hpp \
	inverters/twisted_clover/tmclov_frontends.hpp

#compile or not bgq source
if BGQ
__top_builddir__lib_libnissa_a_SOURCES+= \
	inverters/staggered/cgm_invert_stD2ee_m2_bgq.cpp \
	inverters/staggered/cgm_32_invert_stD2ee_m2_bgq.cpp \
	inverters/staggered/cg_invert_stD2ee_m2_bgq.cpp \
	inverters/twisted_mass/cgm_invert_tmQ2_bgq.cpp \
	inverters/twisted_mass/cg_64_invert_tmQ2_bgq.cpp \
	inverters/twisted_mass/cg_128_invert_tmQ2_bgq.cpp \
	inverters/twisted_mass/cg_64_invert_tmD_eoprec_bgq.cpp \
	inverters/twisted_mass/cg_invert_tmQ2_bgq.cpp \
	inverters/twisted_clover/cgm_invert_tmclovQ2_bgq.cpp \
	inverters/twisted_clover/cg_64_invert_tmclovQ2_bgq.cpp \
	inverters/twisted_clover/cg_64_invert_tmclovD_eoprec_bgq.cpp \
	inverters/twisted_clover/cg_128_invert_tmclovQ2_bgq.cpp

include_HEADERS+= \
	inverters/staggered/cgm_invert_stD2ee_m2_bgq.hpp \
	inverters/staggered/cgm_32_invert_stD2ee_m2_bgq.hpp \
	inverters/staggered/cg_invert_stD2ee_m2_bgq.hpp \
	inverters/twisted_mass/cgm_invert_tmQ2_bgq.hpp \
	inverters/twisted_mass/cg_64_invert_tmQ2_bgq.hpp \
	inverters/twisted_mass/cg_128_invert_tmQ2_bgq.hpp \
	inverters/twisted_mass/cg_64_invert_tmD_eoprec_bgq.hpp \
	inverters/twisted_mass/cg_invert_tmQ2_bgq.hpp \
	inverters/twisted_clover/cg_64_invert_tmclovQ2_bgq.hpp \
	inverters/twisted_clover/cg_64_invert_tmclovD_eoprec_bgq.hpp \
	inverters/twisted_clover/cg_128_invert_tmclovQ2_bgq.hpp \
	inverters/twisted_clover/cgm_invert_tmclovQ2_bgq.hpp \
	inverters/twisted_clover/cgm_invert_tmclovDkern_eoprec_square_portable.hpp
endif

EXTRA_DIST+= \
	inverters/twisted_mass/cg_invert_tmQ2_common.cpp \
	inverters/twisted_mass/cgm_invert_tmQ2_bgq.hpp \
	inverters/templates/cgm_32_invert_template_threaded.cpp \
	inverters/templates/cg_128_invert_template_threaded.cpp \
	inverters/templates/cgm_invert_template_summsol_threaded.hpp \
	inverters/templates/cg_invert_template_threaded.cpp \
	inverters/templates/bicgstab_invert_template_threaded.cpp \
	inverters/templates/cgm_invert_template_threaded.cpp

############################################### io #################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	io/checksum.cpp \
	io/endianness.cpp \
	io/input.cpp \
	io/ILDG_File.cpp \
	io/reader.cpp \
	io/writer.cpp

include_HEADERS+= \
	io/buffer.hpp \
	io/checksum.hpp \
	io/endianness.hpp \
	io/input.hpp \
	io/ILDG_File.hpp \
	io/reader.hpp \
	io/writer.hpp

########################################### linalgs ##################################

__top_builddir__lib_libnissa_a_SOURCES+=linalgs/linalgs.cpp
include_HEADERS+=linalgs/linalgs.hpp

########################################## new types #################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	new_types/complex.cpp \
	new_types/dirac.cpp \
	new_types/dsfmt.cpp \
	new_types/float_128.cpp \
	new_types/float_256.cpp \
	new_types/metadynamics.cpp \
	new_types/high_prec.cpp \
	new_types/rat_approx.cpp \
	new_types/read_new_types.cpp \
	new_types/spin.cpp \
	new_types/su3_op.cpp

include_HEADERS+= \
	new_types/complex.hpp \
	new_types/dirac.hpp \
	new_types/dsfmt.hpp \
	new_types/float_128.hpp \
	new_types/float_256.hpp \
	new_types/high_prec.hpp \
	new_types/metadynamics.hpp \
	new_types/rat_approx.hpp \
	new_types/read_new_types.hpp \
	new_types/spin.hpp \
	new_types/su3.hpp \
	new_types/su3_op.hpp

####################################### operations ##################################

__top_builddir__lib_libnissa_a_SOURCES+= \
	operations/contract/mesons_2pts.cpp \
	operations/contract/mesons_eight.cpp \
	operations/contract/optimized_mesons_2pts.cpp \
	operations/contract/site_contract.cpp \
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
	operations/smearing/gaussian.cpp \
	operations/smearing/HYP.cpp \
	operations/smearing/smooth.cpp \
	operations/smearing/stout_common.cpp \
	operations/smearing/stout_eo.cpp \
	operations/smearing/stout_lx.cpp \
	operations/smearing/Wflow.cpp \
	operations/stag/Laplace_op_2links.cpp \
	operations/stag/mesons.cpp \
	operations/stag/magnetization.cpp \
	operations/stag/nucleon.cpp \
	operations/stag/putpourri.cpp \
	operations/stag/qed_corr.cpp \
	operations/stag/rendens.cpp \
	operations/stag/spinpol.cpp \
	operations/stag/stag.cpp \
	operations/stag/zumba.cpp \
	operations/su3_paths/all_rectangles.cpp \
	operations/su3_paths/arbitrary.cpp \
	operations/su3_paths/clover_term.cpp \
	operations/su3_paths/gauge_sweeper.cpp \
	operations/su3_paths/plaquette.cpp \
	operations/su3_paths/pline.cpp \
	operations/su3_paths/rectangles.cpp \
	operations/su3_paths/rectangular_staples.cpp \
	operations/su3_paths/squared_staples.cpp \
	operations/su3_paths/topological_charge.cpp \
	operations/su3_paths/watusso.cpp \
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
	operations/remez/remez_algorithm.hpp \
	operations/contract/mesons_2pts.hpp \
	operations/contract/mesons_eight.hpp \
	operations/contract/optimized_mesons_2pts.hpp \
	operations/contract/site_contract.hpp \
	operations/smearing/APE.hpp \
	operations/smearing/cooling.hpp \
	operations/smearing/gaussian.hpp \
	operations/smearing/HYP.hpp \
	operations/smearing/smooth.hpp \
	operations/smearing/stout.hpp \
	operations/smearing/Wflow.hpp \
	operations/stag/Laplace_op_2links.hpp \
	operations/stag/mesons.hpp \
	operations/stag/magnetization.hpp \
	operations/stag/nucleon.hpp \
	operations/stag/putpourri.hpp \
	operations/stag/qed_corr.hpp \
	operations/stag/rendens.hpp \
	operations/stag/spinpol.hpp \
	operations/stag/stag.hpp \
	operations/stag/zumba.hpp \
	operations/su3_paths/all_rectangles.hpp \
	operations/su3_paths/arbitrary.hpp \
	operations/su3_paths/clover_term.hpp \
	operations/su3_paths/gauge_sweeper.hpp \
	operations/su3_paths/plaquette.hpp \
	operations/su3_paths/pline.hpp \
	operations/su3_paths/rectangles.hpp \
	operations/su3_paths/rectangular_staples.hpp \
	operations/su3_paths/squared_staples.hpp \
	operations/su3_paths/topological_charge.hpp \
	operations/su3_paths/watusso.hpp \
	operations/vector_gather.hpp

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
