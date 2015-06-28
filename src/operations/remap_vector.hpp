#ifndef _REMAP_VECTOR_HPP
#define _REMAP_VECTOR_HPP

#include "communicate/all_to_all.hpp"

namespace nissa
{
  void remap_lx_vector_to_locd(void *out,void *in,int nbytes,int mu);
  void remap_locd_vector_to_lx(void *out,void *in,int nbytes,int mu);
}

#endif
