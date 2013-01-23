#ifndef _ADDITIONAL_VARIABLES
#define _ADDITIONAL_VARIABLES

#include "tmlQCD_types.h"

#ifdef __cplusplus 

extern "C"
{
  int ****g_ipt;
  tmlQCD_su3 **g_gauge_field; 
  int g_update_gauge_copy;
}

#else

extern int ****g_ipt;
extern tmlQCD_su3 **g_gauge_field; 
extern int g_update_gauge_copy;

#endif

#endif
