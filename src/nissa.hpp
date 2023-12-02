#ifndef _NISSA_HPP
#define _NISSA_HPP

//including config.hpp
#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <base/bench.hpp>
#include <base/close.hpp>
#include <base/debug.hpp>
#ifdef USE_EXTERNAL_SOLVER
# include <base/export_conf_to_external_solver.hpp>
#endif
#include <base/old_field.hpp>
#include <base/init.hpp>
#include <base/lattice.hpp>
#include <base/multiGridParams.hpp>
#include <base/random.hpp>
#include <base/sitmo.hpp>
#ifdef USE_TMLQCD
# include <base/tmLQCD_bridge.hpp>
#endif
#ifdef USE_DDALPHAAMG
# include <base/DDalphaAMG_bridge.hpp>
#endif
#ifdef USE_QUDA
 #include <base/quda_bridge.hpp>
#endif
#include <base/vectors.hpp>
#ifdef USE_CUDA
# include <base/cuda.hpp>
#endif

#include <communicate/allToAll.hpp>
#include <communicate/communicate.hpp>

#include <expr/baseComp.hpp>
#include <expr/comp.hpp>
#include <expr/comps.hpp>
#include <expr/compsMerger.hpp>
#include <expr/compRwCl.hpp>
#include <expr/conj.hpp>
#include <expr/cWiseCombine.hpp>
#include <expr/dagger.hpp>
#include <expr/dynamicTens.hpp>
#include <expr/eoField.hpp>
#include <expr/exponentiator.hpp>
#include <expr/field.hpp>
#include <expr/funcExpr.hpp>
#include <expr/indexComputer.hpp>
#include <expr/kronDelta.hpp>
#include <expr/mergedComps.hpp>
#include <expr/mirroredNode.hpp>
#include <expr/prod.hpp>
#include <expr/scalar.hpp>
#include <expr/shift.hpp>
#include <expr/stackTens.hpp>
#include <expr/trace.hpp>
#include <expr/trigonometry.hpp>
#include <expr/transp.hpp>

#include <free_theory/cg_eoprec_twisted_free_operator.hpp>
#include <free_theory/free_theory_types.hpp>
#include <free_theory/free_theory_types_routines.hpp>
#include <free_theory/tlSym_gauge_propagator.hpp>
#include <free_theory/twisted_free_Dirac_eoprec_operator.hpp>

#include <geometry/geometry_eo.hpp>
#include <geometry/geometry_lx.hpp>
#include <geometry/geometry_mix.hpp>

#include <io/buffer.hpp>
#include <io/checksum.hpp>
#include <io/endianness.hpp>
#include <io/ILDG_File.hpp>
#include <io/input.hpp>
#include <io/reader.hpp>
#include <io/writer.hpp>

#include <linalgs/linalgs.hpp>
#include <linalgs/reduce.hpp>

#include <metaprogramming/concepts.hpp>
#include <metaprogramming/constnessChanger.hpp>
#include <metaprogramming/crtp.hpp>
#include <metaprogramming/extent.hpp>
#include <metaprogramming/feature.hpp>
#include <metaprogramming/hasMember.hpp>
#include <metaprogramming/inline.hpp>

#include <new_types/complex.hpp>
#include <new_types/dirac.hpp>
#include <new_types/float_128.hpp>
#include <new_types/float128class.hpp>
#include <new_types/float_256.hpp>
#include <new_types/high_prec.hpp>
#include <new_types/metadynamics.hpp>
#include <new_types/rat_approx.hpp>
#include <new_types/rng.hpp>
#include <new_types/spin.hpp>
#include <new_types/su3.hpp>

#include <operations/covariant_derivative.hpp>
#include <operations/fft.hpp>
#include <operations/gauge_fixing.hpp>
#include <operations/gaugeconf.hpp>
#include <operations/localizer.hpp>
#include <operations/remap_vector.hpp>
#include <operations/remez/remez_algorithm.hpp>
#include <operations/shift.hpp>
#include <operations/source.hpp>

#include <operations/vector_gather.hpp>

#include <routines/ios.hpp>
#include <routines/math_routines.hpp>
#include <routines/mpi_routines.hpp>

#include <threads/threads.hpp>

#include <tuples/invokeWithTypesOfTuple.hpp>
#include <tuples/tupleDescribe.hpp>
#include <tuples/tupleFilter.hpp>
#include <tuples/typePosition.hpp>

#endif
