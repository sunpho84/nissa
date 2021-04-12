#ifndef _DIRAC_OPERATOR_WSTAT_HPP
#define _DIRAC_OPERATOR_WSTAT_HPP

#include <geometry/geometry_lx.hpp>
#include <new_types/su3.hpp>

namespace nissa
{
  ///Apply the static operator to a spincolor
  void apply_Wstat(spincolor *out,quad_su3 *conf,spincolor *in,const Direction&,int);

  void apply_Wstat2(spincolor *out,quad_su3 *conf,spincolor *temp,spincolor *in);
  void reconstruct_Wstat(spincolor *outminus,spincolor *outplus,quad_su3 *conf,spincolor *in);
}

#endif
