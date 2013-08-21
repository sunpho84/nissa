AUTOMAKE_OPTIONS=subdir-objects

#compile or not spi stuff
if SPI
spi_src_file=\
	src/bgq/spi.cpp \
	src/bgq/spi.h
else
spi_src_file=
endif

#compile or not thread support
if USE_THREADS
thread_src=src/routines/thread.cpp
thread_hdr=src/routines/thread.h
endif

#compile or not virtual parallelization
if USE_VNODES
vnodes_src_file=src/geometry/geometry_vir.h \
	src/geometry/geometry_vir.cpp
endif

#compile or not bgq test and source
if BGQ

#test target
test_bgq_targ=projects/semileptonic/test_bgq
test_bgq_src=projects/semileptonic/test_bgq.cpp

#bgq barrier
if BGQ_EMU
bgq_barrier_src_file=
else
bgq_barrier_src_file=src/bgq/bgq_barrier.c src/bgq/bgq_barrier.h
endif

#source
bgq_src_file=$(bgq_barrier_src_file) \
	src/bgq/Wilson_hopping_matrix_lx_bgq.cpp \
	src/bgq/Wilson_hopping_matrix_lx_bgq.h \
	src/bgq/staggered_hopping_matrix_eo_or_oe_bgq.cpp \
	src/bgq/staggered_hopping_matrix_eo_or_oe_bgq.h \
	src/bgq/bgq_macros.h \
	src/bgq/intrinsic/declare.h \
	src/bgq/intrinsic/load.h \
	src/bgq/intrinsic/mergesplit.h \
	src/bgq/intrinsic/oper.h \
	src/bgq/intrinsic/store.h \
	src/bgq/intrinsic.h \
	src/dirac_operators/dirac_operator_tmQ/dirac_operator_tmQ_bgq.cpp \
	src/dirac_operators/dirac_operator_tmQ/dirac_operator_tmQ_bgq.h \
	src/dirac_operators/dirac_operator_tmQ2/dirac_operator_tmQ2_bgq.cpp \
	src/dirac_operators/dirac_operator_tmQ2/dirac_operator_tmQ2_bgq.h \
	src/inverters/twisted_mass/cgm_invert_tmQ2_bgq.cpp \
	src/inverters/twisted_mass/cgm_invert_tmQ2_bgq.h
else
bgq_src_file=
endif

#add where to find the library and tell to link
AM_LDFLAGS=-L@top_builddir@
LDADD=libnissa.a

#target library nissa
lib_LIBRARIES=libnissa.a

#add macros
AM_CPPFLAGS=@DEFS@

#headers and files for production of nissa
nissa_headers=src/base/close.h \
	src/base/debug.h \
	src/base/global_variables.h \
	src/base/init.h \
	src/base/random.h \
	src/base/vectors.h \
	src/communicate/borders.h \
	src/communicate/communicate.h \
	src/communicate/edges.h \
	src/dirac_operators/dirac_operator_stD/dirac_operator_stD.h \
	src/dirac_operators/dirac_operator_tmDeoimpr/dirac_operator_tmDeoimpr.h \
	src/dirac_operators/dirac_operator_tmQ/dirac_operator_tmQ.h \
	src/dirac_operators/dirac_operator_tmQ_left/dirac_operator_tmQ_left.h \
	src/dirac_operators/dirac_operator_tmQ2/dirac_operator_tmQ2.h \
	src/dirac_operators/dirac_operator_tmQ/reconstruct_tm_doublet.h \
	src/dirac_operators/dirac_operator_tmclovQ/reconstruct_tmclov_doublet.h \
	src/dirac_operators/dirac_operator_Wstat/dirac_operator_Wstat.h \
	src/dirac_operators/dirac_operator_WclovQ/dirac_operator_WclovQ.h \
	src/dirac_operators/dirac_operator_WclovQ2/dirac_operator_WclovQ2.h \
	src/dirac_operators/dirac_operator_tmclovQ/dirac_operator_tmclovQ.h \
	src/dirac_operators/dirac_operator_tmclovQ2/dirac_operator_tmclovQ2.h \
	src/dirac_operators/dirac_operator_tmQ/dirac_operator_tmQ_128.h \
	src/dirac_operators/dirac_operator_tmDeoimpr/dirac_operator_tmDeoimpr_128.h \
	src/dirac_operators/dirac_operator_tmQ2/dirac_operator_tmQ2_128.h \
	src/geometry/geometry_eo.h \
	src/geometry/geometry_lx.h \
	src/geometry/geometry_mix.h \
	src/geometry/geometry_Wsklx.h \
	src/hmc/backfield.h \
	src/hmc/gauge/gluonic_force.h \
	src/hmc/gauge/tree_level_Symanzik_force.h \
	src/hmc/gauge/tree_level_Symanzik_action.h \
	src/hmc/gauge/Wilson_force.h \
	src/hmc/momenta/momenta_generation.h \
	src/hmc/momenta/momenta_action.h \
	src/hmc/rootst_eoimpr/adaptative_control.h \
	src/hmc/rootst_eoimpr/rootst_eoimpr_action.h \
	src/hmc/rootst_eoimpr/rootst_eoimpr_eigenvalues.h \
	src/hmc/rootst_eoimpr/rootst_eoimpr_quark_force.h \
	src/hmc/rootst_eoimpr/rootst_eoimpr_omelyan_integrator.h \
	src/hmc/rootst_eoimpr/rootst_eoimpr_pseudofermions_generation.h \
	src/hmc/rootst_eoimpr/rootst_eoimpr_rhmc_step.h \
	src/inverters/twisted_mass/cg_invert_tmDeoimpr.h \
	src/inverters/twisted_mass/cg_invert_tmQ2.h \
	src/inverters/twisted_mass/cg_64_invert_tmQ2.h \
	src/inverters/twisted_mass/cgm_invert_tmQ2.h \
	src/inverters/twisted_mass/cg_128_invert_tmQ2.h \
	src/inverters/twisted_mass/cg_128_invert_tmQ2.h \
	src/inverters/twisted_mass/cg_128_invert_tmDeoimpr.h \
	src/inverters/twisted_mass/tm_frontends.h \
	src/inverters/staggered/cgm_invert_stD2ee_m2.h \
	src/inverters/staggered/cg_invert_stD2ee_m2.h \
	src/inverters/staggered/cg_invert_stD.h \
	src/inverters/staggered/cg_invert_evn_stD.h \
	src/inverters/Wstat/cg_invert_Wstat.h \
	src/inverters/Wclov/cg_invert_WclovQ.h \
	src/inverters/Wclov/cg_invert_WclovQ2.h \
	src/inverters/tmclov/cg_invert_tmclovQ.h \
	src/inverters/tmclov/cg_invert_tmclovQ2.h \
	src/inverters/tmclov/cgm_invert_tmclovQ2.h \
	src/inverters/tmclov/tmclov_frontends.h \
	src/IO/checksum.h \
	src/IO/endianess.h \
	src/IO/ILDG_File.h \
	src/IO/input.h \
	src/IO/reader.h \
	src/IO/writer.h \
	src/linalgs/linalgs.h \
	src/new_types/complex.h \
	src/new_types/dirac.h \
	src/new_types/float128.h \
	src/new_types/float256.h \
	src/new_types/rat_approx.h \
	src/new_types/read_new_types.h \
	src/new_types/spin.h \
	src/new_types/su3.h \
	src/operations/contract/mesons_2pts.h \
	src/operations/contract/mesons_eight.h \
	src/operations/contract/optimized_mesons_2pts.h \
	src/operations/contract/site_contract.h \
	src/operations/contract/stag.h \
	src/operations/covariant_derivative.h \
	src/operations/fft.h \
	src/operations/fourier_transform.h \
	src/operations/gauge_fixing.h \
	src/operations/gaugeconf.h \
	src/operations/remap_vector.h \
	src/operations/remez/remez_algorithm.h \
	src/operations/shift.h \
	src/operations/source.h \
	src/operations/smearing/APE.h \
	src/operations/smearing/HYP.h \
	src/operations/smearing/stout.h \
	src/operations/smearing/gaussian.h \
	src/operations/su3_paths/all_rectangles.h \
	src/operations/su3_paths/arbitrary.h \
	src/operations/su3_paths/plaquette.h \
	src/operations/su3_paths/pline.h \
	src/operations/su3_paths/rectangles.h \
	src/operations/su3_paths/rectangular_staples.h \
	src/operations/su3_paths/shortest_hypercubic_paths.h \
	src/operations/su3_paths/squared_staples.h \
	src/operations/su3_paths/topological_charge.h \
	src/operations/vector_gather.h \
	src/routines/ios.h \
	src/routines/math.h \
	src/routines/mpi.h \
	$(thread_hdr)
nissa_source_files=src/base/close.cpp \
	src/base/debug.cpp \
	src/base/global_variables.cpp \
	src/base/init.cpp \
	src/base/random.cpp \
	src/base/vectors.cpp \
	src/communicate/borders.cpp \
	src/communicate/edges.cpp \
	src/dirac_operators/dirac_operator_stD/dirac_operator_stD.cpp \
	src/dirac_operators/dirac_operator_tmDeoimpr/dirac_operator_tmDeoimpr.cpp \
	src/dirac_operators/dirac_operator_tmQ/dirac_operator_tmQ.cpp \
	src/dirac_operators/dirac_operator_tmQ_left/dirac_operator_tmQ_left.cpp \
	src/dirac_operators/dirac_operator_tmQ2/dirac_operator_tmQ2.cpp \
	src/dirac_operators/dirac_operator_tmQ/reconstruct_tm_doublet.cpp \
	src/dirac_operators/dirac_operator_tmclovQ/reconstruct_tmclov_doublet.cpp \
	src/dirac_operators/dirac_operator_Wstat/dirac_operator_Wstat.cpp \
	src/dirac_operators/dirac_operator_WclovQ/dirac_operator_WclovQ.cpp \
	src/dirac_operators/dirac_operator_WclovQ2/dirac_operator_WclovQ2.cpp \
	src/dirac_operators/dirac_operator_tmclovQ/dirac_operator_tmclovQ.cpp \
	src/dirac_operators/dirac_operator_tmclovQ2/dirac_operator_tmclovQ2.cpp \
	src/dirac_operators/dirac_operator_tmQ/dirac_operator_tmQ_128.cpp \
	src/dirac_operators/dirac_operator_tmDeoimpr/dirac_operator_tmDeoimpr_128.cpp \
	src/dirac_operators/dirac_operator_tmQ2/dirac_operator_tmQ2_128.cpp \
	src/geometry/geometry_eo.cpp \
	src/geometry/geometry_lx.cpp \
	src/geometry/geometry_mix.cpp \
	src/geometry/geometry_Wsklx.cpp \
	src/hmc/backfield.cpp \
	src/hmc/gauge/gluonic_force.cpp \
	src/hmc/gauge/tree_level_Symanzik_force.cpp \
	src/hmc/gauge/tree_level_Symanzik_action.cpp \
	src/hmc/gauge/Wilson_force.cpp \
	src/hmc/momenta/momenta_generation.cpp \
	src/hmc/momenta/momenta_action.cpp \
	src/hmc/rootst_eoimpr/adaptative_control.cpp \
	src/hmc/rootst_eoimpr/rootst_eoimpr_action.cpp \
	src/hmc/rootst_eoimpr/rootst_eoimpr_eigenvalues.cpp \
	src/hmc/rootst_eoimpr/rootst_eoimpr_quark_force.cpp \
	src/hmc/rootst_eoimpr/rootst_eoimpr_omelyan_integrator.cpp \
	src/hmc/rootst_eoimpr/rootst_eoimpr_pseudofermions_generation.cpp \
	src/hmc/rootst_eoimpr/rootst_eoimpr_rhmc_step.cpp \
	src/inverters/twisted_mass/cg_invert_tmDeoimpr.cpp \
	src/inverters/twisted_mass/cg_invert_tmQ2.cpp \
	src/inverters/twisted_mass/cg_64_invert_tmQ2.cpp \
	src/inverters/twisted_mass/cgm_invert_tmQ2.cpp \
	src/inverters/twisted_mass/cg_128_invert_tmQ2.cpp \
	src/inverters/twisted_mass/cg_128_invert_tmQ2.cpp \
	src/inverters/twisted_mass/cg_128_invert_tmDeoimpr.cpp \
	src/inverters/twisted_mass/tm_frontends.cpp \
	src/inverters/staggered/cgm_invert_stD2ee_m2.cpp \
	src/inverters/staggered/cg_invert_stD2ee_m2.cpp \
	src/inverters/staggered/cg_invert_stD.cpp \
	src/inverters/staggered/cg_invert_evn_stD.cpp \
	src/inverters/Wstat/cg_invert_Wstat.cpp \
	src/inverters/Wclov/cg_invert_WclovQ.cpp \
	src/inverters/Wclov/cg_invert_WclovQ2.cpp \
	src/inverters/tmclov/cg_invert_tmclovQ.cpp \
	src/inverters/tmclov/cg_invert_tmclovQ2.cpp \
	src/inverters/tmclov/cgm_invert_tmclovQ2.cpp \
	src/inverters/tmclov/tmclov_frontends.cpp \
	src/IO/checksum.cpp \
	src/IO/endianess.cpp \
	src/IO/ILDG_File.cpp \
	src/IO/input.cpp \
	src/IO/reader.cpp \
	src/IO/writer.cpp \
	src/linalgs/linalgs.cpp \
	src/new_types/complex.cpp \
	src/new_types/dirac.cpp \
	src/new_types/float128.cpp \
	src/new_types/float256.cpp \
	src/new_types/rat_approx.cpp \
	src/new_types/read_new_types.cpp \
	src/new_types/spin.cpp \
	src/new_types/su3.cpp \
	src/operations/contract/mesons_2pts.cpp \
	src/operations/contract/mesons_eight.cpp \
	src/operations/contract/optimized_mesons_2pts.cpp \
	src/operations/contract/site_contract.cpp \
	src/operations/contract/stag.cpp \
	src/operations/covariant_derivative.cpp \
	src/operations/fft.cpp \
	src/operations/fourier_transform.cpp \
	src/operations/gauge_fixing.cpp \
	src/operations/gaugeconf.cpp \
	src/operations/remap_vector.cpp \
	src/operations/remez/remez_algorithm.cpp \
	src/operations/shift.cpp \
	src/operations/source.cpp \
	src/operations/smearing/APE.cpp \
	src/operations/smearing/HYP.cpp \
	src/operations/smearing/stout.cpp \
	src/operations/smearing/gaussian.cpp \
	src/operations/su3_paths/all_rectangles.cpp \
	src/operations/su3_paths/arbitrary.cpp \
	src/operations/su3_paths/plaquette.cpp \
	src/operations/su3_paths/pline.cpp \
	src/operations/su3_paths/rectangles.cpp \
	src/operations/su3_paths/rectangular_staples.cpp \
	src/operations/su3_paths/shortest_hypercubic_paths.cpp \
	src/operations/su3_paths/squared_staples.cpp \
	src/operations/su3_paths/topological_charge.cpp \
	src/operations/vector_gather.cpp \
	src/routines/ios.cpp \
	src/routines/math.cpp \
	src/routines/mpi.cpp \
	$(thread_src)
extra_headers=src/new_types/new_types_definitions.h src/base/macros.h src/base/thread_macros.h

#instruct how to build the library
libnissa_a_SOURCES=$(nissa_source_files) $(nissa_headers) $(extra_headers) \
	$(vnodes_src_file) $(bgq_src_file) $(exp_bgq_src_file) $(spi_src_file)

#########################################################################################################################

#list of all the programs
bin_PROGRAMS=\
	projects/generate_confs/generate_stag \
	$(test_bgq_targ) \
	projects/semileptonic/test_Wsk \
	projects/semileptonic/test_cov \
	projects/semileptonic/test_proj \
	projects/semileptonic/semileptonic_smeared \
	projects/semileptonic/semileptonic_smeared_point_source \
	projects/reno_const/RIMOM \
	projects/eight_BK/smeared_BK_all_in_one \
	tools/unitarity_check/unitarity_check

#	projects/semileptonic/test 
#	projects/static_potential/compute_potential

#program to generate confs
projects_generate_confs_generate_stag_SOURCES=projects/generate_confs/generate_stag.cpp $(extra_headers)

#program to compute BK
projects_eight_BK_smeared_BK_all_in_one_SOURCES=projects/eight_BK/smeared_BK_all_in_one.cpp $(extra_headers)

#program for semileptonic smeared 2 and 3 points
projects_semileptonic_semileptonic_smeared_SOURCES=projects/semileptonic/semileptonic_smeared.cpp $(extra_headers)
projects_semileptonic_test_cov_SOURCES=projects/semileptonic/test_cov.cpp $(extra_headers)
projects_semileptonic_test_Wsk_SOURCES=projects/semileptonic/test_Wsk.cpp $(extra_headers)
if BGQ
projects_semileptonic_test_bgq_SOURCES=projects/semileptonic/test_bgq.cpp $(extra_headers)
endif
projects_semileptonic_test_proj_SOURCES=projects/semileptonic/test_proj.cpp $(extra_headers)
projects_semileptonic_semileptonic_smeared_point_source_SOURCES=projects/semileptonic/semileptonic_smeared_point_source.cpp $(extra_headers)
#projects_semileptonic_test_SOURCES=projects/semileptonic/test.cpp projects/semileptonic/new_types.cpp projects/semileptonic/new_types.h $(extra_headers)
#projects_semileptonic_semileptonic_smeared_point_source_SOURCES=projects/semileptonic/semileptonic_smeared_point_source.cpp $(extra_headers)

#compute static potential
#projects_static_potential_compute_potential_SOURCES=projects/static_potential/compute_potential.cpp $(extra_headers)

#unitarity check
tools_unitarity_check_unitarity_check_SOURCES=tools/unitarity_check/unitarity_check.cpp $(extra_headers)

#RIMOM
projects_reno_const_RIMOM_SOURCES=projects/reno_const/RIMOM.cpp $(extra_headers)

#work_in_prog_paths_su3_paths_SOURCES=work_in_prog/paths/su3_paths.cpp $(extra_headers)

###########################################################################################################################################################

#adds extra files to distribution
EXTRA_DIST= \
	svn_revision \
	README \
	src/dirac_operators/dirac_operator_WclovQ/dirac_operator_WclovQ_portable.cpp \
	src/dirac_operators/dirac_operator_stD/dirac_operator_stD_portable.cpp \
	src/dirac_operators/dirac_operator_tmDeoimpr/dirac_operator_tmDeoimpr_portable.cpp \
	src/dirac_operators/dirac_operator_tmQ/dirac_operator_tmQ_portable.cpp \
	src/dirac_operators/dirac_operator_tmQ_left/dirac_operator_tmQ_left_portable.cpp \
	src/hmc/rootst_eoimpr/rat_expansion_database.cpp \
	src/inverters/twisted_mass/cg_invert_tmQ2_common.cpp \
	src/nissa.h \
	src/create_headers.sh \
	src/inverters/templates/cg_128_invert_template.cpp \
	src/inverters/templates/cg_invert_template.cpp \
	src/inverters/templates/cgm_invert_template.cpp \
	src/inverters/templates/cg_128_invert_template_threaded.cpp \
	src/inverters/templates/cg_invert_template_threaded.cpp \
	src/inverters/templates/cgm_invert_template_threaded.cpp \
	tools/clusterize.cpp \
	tools/clusterize2.cpp \
	tools/pack.c


