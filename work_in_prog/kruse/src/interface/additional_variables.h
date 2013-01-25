#ifndef _ADDITIONAL_VARIABLES
#define _ADDITIONAL_VARIABLES

#include "tmlQCD_types.h"

#ifdef __cplusplus 

extern "C"
{
  int ****g_ipt;
  tmlQCD_su3 **g_gauge_field; 
  int g_update_gauge_copy;
  //tobeadded
  _Complex double ka0,ka1,ka2,ka3;
  int *g_eo2lexic,*g_lexic2eo,*g_lexic2eosub;
}

#else

extern int ****g_ipt;
extern tmlQCD_su3 **g_gauge_field; 
extern int g_update_gauge_copy;
extern _Complex double ka0,ka1,ka2,ka3;
extern int *g_eo2lexic,*g_lexic2eo,*g_lexic2eosub;

#endif

#endif
