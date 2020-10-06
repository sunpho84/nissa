#ifndef _NISSA_HPP
#define _NISSA_HPP

//including config.hpp
#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/bench.hpp"
#include "base/close.hpp"
#include "base/debug.hpp"
#include "base/init.hpp"
#include "base/random.hpp"
#include "base/thread_macros.hpp"
#ifdef USE_TMLQCD
 #include "base/tmLQCD_bridge.hpp"
#endif
#ifdef USE_DDALPHAAMG
 #include "base/DDalphaAMG_bridge.hpp"
#endif
#include "base/vectors.hpp"

//include bg/q specifics
#ifdef BGQ
 #include "bgq/intrinsic.hpp"
 #include "dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_bgq.hpp"
 #include "dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec_bgq.hpp"
 #include "dirac_operators/tmQ/dirac_operator_tmQ_bgq.hpp"
 #include "dirac_operators/tmQ2/dirac_operator_tmQ2_bgq.hpp"
 #include "dirac_operators/stD/dirac_operator_stD_bgq.hpp"
#endif
#ifdef SPI
 #include "bgq/spi.hpp"
#endif

#ifdef USE_VNODES
 #include "geometry/geometry_vir.hpp"
#endif

#include "communicate/borders.hpp"
#include "communicate/communicate.hpp"
#include "communicate/edges.hpp"

#include "dirac_operators/momenta/MFACC.hpp"
#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec.hpp"
#include "dirac_operators/tmQ/dirac_operator_tmQ.hpp"
#include "dirac_operators/tmQ/dirac_operator_tmQ_128.hpp"
#include "dirac_operators/tmQ/reconstruct_tm_doublet.hpp"
#include "dirac_operators/tmQ2/dirac_operator_tmQ2.hpp"
#include "dirac_operators/tmQ2/dirac_operator_tmQ2_128.hpp"
#include "dirac_operators/tmQ_left/dirac_operator_tmQ_left.hpp"
#include "dirac_operators/Wstat/dirac_operator_Wstat.hpp"
#include "dirac_operators/WclovQ/dirac_operator_WclovQ.hpp"
#include "dirac_operators/WclovQ2/dirac_operator_WclovQ2.hpp"
#include "dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec.hpp"
#include "dirac_operators/tmclovQ/dirac_operator_tmclovQ.hpp"
#include "dirac_operators/tmclovQ/reconstruct_tmclov_doublet.hpp"
#include "dirac_operators/tmclovQ2/dirac_operator_tmclovQ2.hpp"

#include "free_theory/cg_eoprec_twisted_free_operator.hpp"
#include "free_theory/free_theory_types.hpp"
#include "free_theory/free_theory_types_routines.hpp"
#include "free_theory/tlSym_gauge_propagator.hpp"
#include "free_theory/twisted_propagator.hpp"
#include "free_theory/twisted_free_Dirac_eoprec_operator.hpp"

#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_Leb.hpp"
#include "geometry/geometry_mix.hpp"

#include "hmc/backfield.hpp"
#include "hmc/fermions/rootst_eoimpr_quark_force.hpp"
#include "hmc/fermions/pseudofermions_generation.hpp"
#include "hmc/gauge/gluonic_action.hpp"
#include "hmc/gauge/gluonic_force.hpp"
#include "hmc/gauge/pure_gauge_hmc_step.hpp"
#include "hmc/gauge/pure_gauge_Omelyan_integrator.hpp"
#include "hmc/gauge/pure_gauge_implicit_integrator.hpp"
#include "hmc/gauge/topological_action.hpp"
#include "hmc/gauge/topological_force.hpp"
#include "hmc/gauge/Symanzik_force.hpp"
#include "hmc/gauge/Symanzik_action.hpp"
#include "hmc/gauge/Wilson_action.hpp"
#include "hmc/gauge/Wilson_force.hpp"
#include "hmc/hmc.hpp"
#include "hmc/momenta/momenta_action.hpp"
#include "hmc/momenta/momenta_evolve.hpp"
#include "hmc/momenta/momenta_generation.hpp"
#include "hmc/multipseudo/Omelyan_integrator.hpp"
#include "hmc/multipseudo/set_expansions.hpp"
#include "hmc/multipseudo/theory_action.hpp"
#include "hmc/multipseudo/multipseudo_rhmc_step.hpp"

#include "inverters/staggered/cg_invert_stD.hpp"
#include "inverters/twisted_mass/cg_128_invert_tmQ2.hpp"
#include "inverters/twisted_mass/cg_invert_tmD_eoprec.hpp"
#include "inverters/twisted_mass/cg_invert_tmQ2.hpp"
#include "inverters/twisted_mass/cgm_invert_tmQ2.hpp"
#include "inverters/twisted_mass/tm_frontends.hpp"
#include "inverters/Wclov/cg_invert_WclovQ2.hpp"
#include "inverters/Wclov/cg_invert_WclovQ.hpp"
#include "inverters/twisted_clover//cg_invert_tmclovD_eoprec.hpp"
#include "inverters/twisted_clover/cg_invert_tmclovQ2.hpp"
#include "inverters/twisted_clover/cg_invert_tmclovQ.hpp"
#include "inverters/twisted_clover/cgm_invert_tmclovQ2.hpp"
#include "inverters/twisted_clover/tmclov_frontends.hpp"

#include "io/buffer.hpp"
#include "io/checksum.hpp"
#include "io/endianness.hpp"
#include "io/ILDG_File.hpp"
#include "io/input.hpp"
#include "io/reader.hpp"
#include "io/writer.hpp"

#include "linalgs/linalgs.hpp"

#include "new_types/complex.hpp"
#include "new_types/dirac.hpp"
#include "new_types/dsfmt.hpp"
#include "new_types/float_128.hpp"
#include "new_types/float_256.hpp"
#include "new_types/high_prec.hpp"
#include "new_types/metadynamics.hpp"
#include "new_types/rat_approx.hpp"
#include "new_types/read_new_types.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3.hpp"

#include "operations/contract/mesons_2pts.hpp"
#include "operations/contract/mesons_eight.hpp"
#include "operations/contract/optimized_mesons_2pts.hpp"
#include "operations/contract/site_contract.hpp"
#include "operations/covariant_derivative.hpp"
#include "operations/fft.hpp"
#include "operations/fourier_transform.hpp"
#include "operations/gauge_fixing.hpp"
#include "operations/gaugeconf.hpp"
#include "operations/remap_vector.hpp"
#include "operations/remez/remez_algorithm.hpp"
#include "operations/shift.hpp"
#include "operations/smearing/APE.hpp"
#include "operations/smearing/HYP.hpp"
#include "operations/smearing/gaussian.hpp"
#include "operations/smearing/recursive_Wflower.hpp"
#include "operations/smearing/stout.hpp"
#include "operations/smearing/smooth.hpp"
#include "operations/smearing/Wflow.hpp"
#include "operations/source.hpp"

#include "operations/stag/magnetization.hpp"
#include "operations/stag/mesons.hpp"
#include "operations/stag/nucleon.hpp"
#include "operations/stag/putpourri.hpp"
#include "operations/stag/qed_corr.hpp"
#include "operations/stag/rendens.hpp"
#include "operations/stag/spinpol.hpp"
#include "operations/stag/stag.hpp"
#include "operations/stag/zumba.hpp"

#include "operations/su3_paths/all_rectangles.hpp"
#include "operations/su3_paths/arbitrary.hpp"
#include "operations/su3_paths/clover_term.hpp"
#include "operations/su3_paths/gauge_sweeper.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "operations/su3_paths/pline.hpp"
#include "operations/su3_paths/rectangles.hpp"
#include "operations/su3_paths/rectangular_staples.hpp"
#include "operations/su3_paths/squared_staples.hpp"
#include "operations/su3_paths/topological_charge.hpp"
#include "operations/su3_paths/watusso.hpp"

#include "operations/vector_gather.hpp"

#include "routines/ios.hpp"
#include "routines/math_routines.hpp"
#include "routines/mpi_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#endif
