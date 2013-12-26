#ifndef _REMAP_VECTOR_H
#define _REMAP_VECTOR_H

#include "communicate/all_to_all.hpp"

namespace nissa
{
  struct vector_remap_t : all_to_all_comm_t
  {
    vector_remap_t(int nel_out,void (*index)(int &irank_to,int &iel_to,int iel_fr,void *pars),void *pars);
    void remap(void *out,void *in,int bps){communicate(out,in,bps);}
  };
}

#endif
