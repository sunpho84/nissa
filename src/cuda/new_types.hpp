#ifndef _CUDA_NEW_TYPES_H
#define _CUDA_NEW_TYPES_H

#include "macros.hpp"

extern DEVICE_CONSTANT int loc_volh;

namespace cuda
{
  template <class T,const int nreals_per_site> class cuda_field
  {
  private:
    T *data;
  public:
    HOST cuda_field() {data=NULL;}
    DEVICE T*& operator[](int id){return data+id*loc_volh;};
    HOST void reset() 
    {
      dim3 block_dimension(NCUDA_THREADS);
      dim3 grid_dimension(nissa::loc_volh/block_dimension.x);

      field_reset<<<grid_dimension,block_dimension>>>(data,nreals_per_site);
    }
  };
  
  class float_gauge_field : cuda_field<float,12>
  {
  };
}

#endif
