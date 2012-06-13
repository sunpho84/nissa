#include "../../../../src/nissa.h"

#include "../types/types.h"

//prop1 will be reverted
void compute_all_2pts_qdagq_correlations(corr16 *corr,spinspin *prop1,spinspin *prop2)
{
  dirac_matr t1[16],t2[16];
  for(int igamma=0;igamma<16;igamma++)
    {
      dirac_prod(&(t1[igamma]),&(base_gamma[igamma]),&(base_gamma[5]));
      dirac_prod(&(t2[igamma]),&(base_gamma[5]),&(base_gamma[igamma]));
    }
  
  nissa_loc_vol_loop(ivol)
    for(int igamma=0;igamma<16;igamma++)
      site_trace_g_sdag_g_s(corr[ivol][igamma],&(t1[igamma]),prop1[ivol],&(t2[igamma]),prop2[ivol]);
}
