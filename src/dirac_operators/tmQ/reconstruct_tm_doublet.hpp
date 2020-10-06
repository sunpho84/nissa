#ifndef _RECONSTRUCT_TM_DOUBLET_HPP
#define _RECONSTRUCT_TM_DOUBLET_HPP

#include "new_types/su3.hpp"


namespace nissa
{
  void reconstruct_tm_doublet(spincolor *outminus,spincolor *outplus,quad_su3 *conf,double kappac,double mu,spincolor *in);
}

#endif
