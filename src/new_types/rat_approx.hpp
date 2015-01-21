#ifndef _RAT_APPROX_HPP
#define _RAT_APPROX_HPP

#include "io/buffer.hpp"

namespace nissa
{
  void convert_rat_approx(buffer_t &s,rat_approx_t *appr,int nflav);
  void convert_rat_approx(char *&data,int &data_length,rat_approx_t *appr,int nflav);
  void convert_rat_approx(rat_approx_t *&appr,int &nflav,char *data,int data_length);
  void master_printf_rat_approx(rat_approx_t *appr);
  void rat_approx_create(rat_approx_t *appr,int nterms,const char *name=NULL);
  void rat_approx_destroy(rat_approx_t *appr);
  //shift all poles
  inline void rat_approx_shift_all_poles(rat_approx_t *appr,double sh)
  {for(int iterm=0;iterm<appr->degree;iterm++) appr->poles[iterm]+=sh;}
  
}

#endif
