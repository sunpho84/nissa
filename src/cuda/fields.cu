#include "macros.hpp"

namespace cuda
{
  template <class T> GLOBAL void field_reset(T *d,const int nper_site) 
  {
    //unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    //const unsigned int grid_length = blockDim.x * gridDim.x;

    //__shared__ double norm[128];  //Allocates shared mem

    //float4 f_0;  
  }
  
  GLOBAL void field_reset(float *d,const int nper_site);
}
