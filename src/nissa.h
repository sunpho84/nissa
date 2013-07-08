#ifndef _NISSA_H
#define _NISSA_H

//including config.h
#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "base/close.h"
#include "base/debug.h"
#include "base/global_variables.h"
#include "base/init.h"
#include "base/macros.h"
#include "base/random.h"
#include "base/thread_macros.h"
#include "base/vectors.h"

//include bg/q specifics
#ifdef BGQ
 #include "bgq/geometry_bgq.h"
 #include "bgq/intrinsic.h"
 #include "dirac_operators/dirac_operator_tmQ/dirac_operator_tmQ_bgq.h"
 #include "dirac_operators/dirac_operator_tmQ2/dirac_operator_tmQ2_bgq.h"
#endif
#ifdef SPI
 #include "bgq/spi.h"
#endif

#include "communicate/communicate.h"

#include "hmc/backfield.h"
#include "hmc/gauge/gluonic_force.h"
#include "hmc/gauge/Wilson_force.h"
#include "hmc/gauge/tree_level_Symanzik_force.h"
#include "hmc/gauge/tree_level_Symanzik_action.h"
#include "hmc/momenta/momenta_action.h"
#include "hmc/momenta/momenta_generation.h"
#include "hmc/rootst_eoimpr/adaptative_control.h"
#include "hmc/rootst_eoimpr/rootst_eoimpr_action.h"
#include "hmc/rootst_eoimpr/rootst_eoimpr_eigenvalues.h"
#include "hmc/rootst_eoimpr/rootst_eoimpr_quark_force.h"
#include "hmc/rootst_eoimpr/rootst_eoimpr_omelyan_integrator.h"
#include "hmc/rootst_eoimpr/rootst_eoimpr_pseudofermions_generation.h"
#include "hmc/rootst_eoimpr/rootst_eoimpr_rhmc_step.h"

#include "IO/checksum.h"
#include "IO/endianess.h"
#include "IO/ILDG_File.h"
#include "IO/input.h"
#include "IO/reader.h"
#include "IO/writer.h"

#include "dirac_operators/dirac_operator_stD/dirac_operator_stD.h"
#include "dirac_operators/dirac_operator_tmDeoimpr/dirac_operator_tmDeoimpr.h"
#include "dirac_operators/dirac_operator_tmQ/dirac_operator_tmQ.h"
#include "dirac_operators/dirac_operator_tmQ/dirac_operator_tmQ_128.h"
#include "dirac_operators/dirac_operator_tmQ/reconstruct_tm_doublet.h"
#include "dirac_operators/dirac_operator_tmQ2/dirac_operator_tmQ2.h"
#include "dirac_operators/dirac_operator_tmQ2/dirac_operator_tmQ2_128.h"
#include "dirac_operators/dirac_operator_tmQ_left/dirac_operator_tmQ_left.h"
#include "dirac_operators/dirac_operator_Wstat/dirac_operator_Wstat.h"
#include "dirac_operators/dirac_operator_WclovQ/dirac_operator_WclovQ.h"
#include "dirac_operators/dirac_operator_WclovQ2/dirac_operator_WclovQ2.h"
#include "dirac_operators/dirac_operator_tmclovQ/dirac_operator_tmclovQ.h"
#include "dirac_operators/dirac_operator_tmclovQ/reconstruct_tmclov_doublet.h"
#include "dirac_operators/dirac_operator_tmclovQ2/dirac_operator_tmclovQ2.h"

#include "geometry/geometry_eo.h"
#include "geometry/geometry_lx.h"
#include "geometry/geometry_mix.h"

#include "inverters/twisted_mass/cg_128_invert_tmQ2.h"
#include "inverters/twisted_mass/cg_invert_tmDeoimpr.h"
#include "inverters/twisted_mass/cg_invert_tmQ2.h"
#include "inverters/twisted_mass/cgm_invert_tmQ2.h"
#include "inverters/twisted_mass/tm_frontends.h"
#include "inverters/Wstat/cg_invert_Wstat.h"
#include "inverters/Wclov/cg_invert_WclovQ2.h"
#include "inverters/Wclov/cg_invert_WclovQ.h"
#include "inverters/tmclov/cg_invert_tmclovQ2.h"
#include "inverters/tmclov/cg_invert_tmclovQ.h"
#include "inverters/tmclov/cgm_invert_tmclovQ2.h"
#include "inverters/tmclov/tmclov_frontends.h"

#include "linalgs/linalgs.h"

#include "new_types/complex.h"
#include "new_types/dirac.h"
#include "new_types/float128.h"
#include "new_types/float256.h"
#include "new_types/new_types_definitions.h"
#include "new_types/rat_approx.h"
#include "new_types/read_new_types.h"
#include "new_types/spin.h"
#include "new_types/su3.h"

#include "operations/contract/mesons_2pts.h"
#include "operations/contract/mesons_eight.h"
#include "operations/contract/site_contract.h"
#include "operations/contract/stag.h"
#include "operations/covariant_derivative.h"
#include "operations/fft.h"
#include "operations/fourier_transform.h"
#include "operations/gauge_fixing.h"
#include "operations/gaugeconf.h"
#include "operations/remap_vector.h"
#include "operations/remez/remez_algorithm.h"
#include "operations/shift.h"
#include "operations/smearing/APE.h"
#include "operations/smearing/HYP.h"
#include "operations/smearing/stout.h"
#include "operations/smearing/gaussian.h"
#include "operations/source.h"

#include "operations/su3_paths/all_rectangles.h"
#include "operations/su3_paths/arbitrary.h"
#include "operations/su3_paths/plaquette.h"
#include "operations/su3_paths/pline.h"
#include "operations/su3_paths/rectangles.h"
#include "operations/su3_paths/squared_staples.h"
#include "operations/su3_paths/shortest_hypercubic_paths.h"
#include "operations/su3_paths/topological_charge.h"

#include "operations/vector_gather.h"

#include "routines/ios.h"
#include "routines/math.h"
#include "routines/mpi.h"
#ifdef USE_THREADS
 #include "routines/thread.h"
#endif

#endif
