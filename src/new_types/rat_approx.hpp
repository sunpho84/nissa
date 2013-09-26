#ifndef _RAT_EXP_H
#define _RAT_EXP_H

namespace nissa
{
  void master_printf_rat_approx(rat_approx_t *appr);
  void rat_approx_create(rat_approx_t *appr,int nterms,const char *name);
  void rat_approx_destroy(rat_approx_t *appr);
}

#endif
