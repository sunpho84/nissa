#ifndef _BARYON_HPP
#define _BARYON_HPP

#include "new_types/su3.hpp"

namespace nissa
{
  /// Compute the barionic contractions, for a specific projection
  ///
  /// The two gammas refers to the source and sink projection. The two
  /// Wick contractions are stored separately. The three quarks are
  /// intended to be l-d-l' (like-dislike-like), and the two Wick
  /// contractions are the direct and exchange respectively. The
  /// second Wick contraction makes sense only if l' and l are the
  /// same. Propagators are expected to contain a phase in time, and
  /// are automatically shifted to zero time. Each propagator is a
  /// list of 12 spincolor vector, each corresponding to the twelve
  /// components.
  void compute_baryon_2pts_proj_contr(complex* contr,             ///< Output, with indices {t,iWick}
				      const int& igSo,            ///< Gamma in the source
				      const int& igSi,            ///< Gamma in the sink
				      spincolor** Q1,             ///< l propagator
				      spincolor** Q2,             ///< d propagator
				      spincolor** Q3,             ///< l' propagator
				      const int source_coord,     ///< Source coordinate
				      const double& temporal_bc); ///< Boundary conditon in time
  
  /// Compute the nucleon contractions
  ///
  /// To compute the proton, pass U as l and D as d. For neutron, the
  /// opposite. The two Wick contractions are internally combined. The
  /// propagators are expected to follow the same convention of above.
  void compute_nucleon_2pts_contr(complex* contr,             ///< Output, with indices {t,iWick}
				  spincolor** Ql,             ///< l propagator
				  spincolor** Qd,             ///< d propagator
				  const int source_coord,     ///< Source coordinate
				  const double& temporal_bc); ///< Boundary conditon in time
}

#endif
