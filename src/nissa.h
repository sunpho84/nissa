#ifndef _NISSA_H
#define _NISSA_H

#include "base/close.h"
#include "base/communicate.h"
#include "base/debug.h"
#include "base/global_variables.h"
#include "base/init.h"
#include "base/macros.h"
#include "base/random.h"
#include "base/routines.h"
#include "base/sse_instructions.h"
#include "base/vectors.h"

#ifdef OMP
 #include <omp.h>
#endif

#ifdef BGP
 #ifndef BGP_EMU
  #include <builtins.h>
 #endif
 #include "base/bgp_instructions.h"
#endif

#include "hmc/backfield.h"
#include "hmc/gauge/wilson_force.h"
#include "hmc/momenta/momenta_action.h"
#include "hmc/momenta/momenta_generation.h"
#include "hmc/rootst_eoimpr/rootst_eoimpr_action.h"
#include "hmc/rootst_eoimpr/rootst_eoimpr_eigenvalues.h"
#include "hmc/rootst_eoimpr/rootst_eoimpr_force.h"
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

#include "geometry/geometry_eo.h"
#include "geometry/geometry_lx.h"
#include "geometry/geometry_mix.h"

#include "inverters/twisted_mass/cg_128_invert_tmQ2.h"
#include "inverters/twisted_mass/cg_invert_tmDeoimpr.h"
#include "inverters/twisted_mass/cg_invert_tmQ2.h"
#include "inverters/twisted_mass/cgm_invert_tmQ2.h"
#include "inverters/twisted_mass/tm_frontends.h"
#include "inverters/Wstat/cg_invert_Wstat.h"

#include "linalgs/linalgs.h"

#include "new_types/complex.h"
#include "new_types/dirac.h"
#include "new_types/float128.h"
#include "new_types/new_types_definitions.h"
#include "new_types/rat_exp.h"
#include "new_types/spin.h"
#include "new_types/su3.h"

#include "operations/contract.h"
#include "operations/contract_stag.h"
#include "operations/fft.h"
#include "operations/fourier_transform.h"
#include "operations/gauge_fixing.h"
#include "operations/gaugeconf.h"
#include "operations/remap_vector.h"
#include "operations/smear.h"
#include "operations/source.h"
#include "operations/su3_paths.h"
#include "operations/vector_gather.h"

#endif
