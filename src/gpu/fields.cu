#include "macros.hpp"

namespace cuda
{
  template <class T> GLOBAL void field_reset(T *d,const int nper_site) 
  {
    unsigned int idx=blockIdx.x*blockDim.x+threadIdx.x;
    for(int i=0;i<nper_site;i++) d[i][idx]=0;
  }
  
  GLOBAL void field_reset(float *d,const int nper_site);
}
