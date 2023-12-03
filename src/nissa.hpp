#ifndef _NISSA_HPP
#define _NISSA_HPP

//including config.hpp
#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <base/close.hpp>
#include <base/debug.hpp>
#ifdef USE_EXTERNAL_SOLVER
# include <base/export_conf_to_external_solver.hpp>
#endif
#include <base/init.hpp>
#include <base/lattice.hpp>
#include <base/multiGridParams.hpp>
#ifdef USE_DDALPHAAMG
# include <base/DDalphaAMG_bridge.hpp>
#endif
#ifdef USE_QUDA
 #include <base/quda_bridge.hpp>
#endif
#include <base/qcd.hpp>
#include <base/vectors.hpp>
#ifdef USE_CUDA
# include <base/cuda.hpp>
#endif

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

#include <geometry/geometry_eo.hpp>
#include <geometry/geometry_lx.hpp>

#include <io/buffer.hpp>
#include <io/checksum.hpp>
#include <io/endianness.hpp>
#include <io/ILDG_File.hpp>
#include <io/reader.hpp>
#include <io/writer.hpp>

#include <linalgs/reduce.hpp>

#include <metaprogramming/concepts.hpp>
#include <metaprogramming/constnessChanger.hpp>
#include <metaprogramming/crtp.hpp>
#include <metaprogramming/extent.hpp>
#include <metaprogramming/feature.hpp>
#include <metaprogramming/hasMember.hpp>
#include <metaprogramming/inline.hpp>

#include <new_types/float128class.hpp>
#include <new_types/rng.hpp>

#include <operations/fft.hpp>
#include <operations/localizer.hpp>

#include <routines/ios.hpp>
#include <routines/math_routines.hpp>
#include <routines/mpi_routines.hpp>

#include <threads/threads.hpp>

#include <tuples/invokeWithTypesOfTuple.hpp>
#include <tuples/tupleDescribe.hpp>
#include <tuples/tupleFilter.hpp>
#include <tuples/typePosition.hpp>

#endif
