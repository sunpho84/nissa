#ifndef _FREE_THEORY_TYPES_HPP
#define _FREE_THEORY_TYPES_HPP

#include <stdint.h>

#include "geometry/geometry_lx.hpp"
#include "new_types/complex.hpp"

namespace nissa
{
  enum tm_basis_t{WILSON_BASE,MAX_TWIST_BASE};
  enum zero_mode_sub_t{PECIONA,UNNO_ALEMANNA,ONLY_100};
  const double FEYNMAN_ALPHA=1,LANDAU_ALPHA=0;
  const double WILSON_C1=0,TLSYM_C1=-1.0/12;
  
  /// Information on a twisted mass quark
  struct tm_quark_info
  {
    double kappa;
    double mass;
    Momentum bc;
    double zmp;
    int r;
    
    // Once tensor has a copy constructor, remove these
    
    /// Copy constructor
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    tm_quark_info(const tm_quark_info& oth)
    {
      kappa=oth.kappa;
      mass=oth.mass;
      
      FOR_ALL_DIRS(mu)
	bc(mu)=oth.bc(mu);
      
      zmp=oth.zmp;
      r=oth.r;
    }
    
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    tm_quark_info(const double& kappa,const double& mass,const int& r,const double& theta) :
      kappa(kappa),mass(mass),zmp(0),r(r)
    {
      bc(tDir)=1;
      for(Dir mu=1;mu<NDIM;mu++)
	bc(mu)=theta;
    }
    
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    tm_quark_info(const double& kappa,const double& mass,const int& r,const Momentum& _bc) :
      kappa(kappa),mass(mass),zmp(0),r(r)
    {
      FOR_ALL_DIRS(mu)
	bc(mu)=_bc(mu);
    }
    
    tm_quark_info() {}
  };
  
  struct gauge_info
  {
    zero_mode_sub_t zms;
    double alpha;
    double c1;
    Momentum bc;
    
    gauge_info()
    {
      zms=UNNO_ALEMANNA;
      alpha=LANDAU_ALPHA;
      c1=WILSON_C1;
      
      for(Dir mu=0;mu<NDIM;mu++)
	bc(mu)=0;
    }
    
    CUDA_HOST_AND_DEVICE gauge_info(const gauge_info& oth) //nasty
    {
      zms=oth.zms;
      alpha=oth.alpha;
      c1=oth.c1;
      
      FOR_ALL_DIRS(mu)
	bc(mu)=oth.bc(mu);
    }
};
  
  typedef complex corr16[16];
  
  struct bmpfile_magic{unsigned char magic[2];};
  
  struct bmpfile_header
  {
    uint32_t filesz;
    uint16_t creator1;
    uint16_t creator2;
    uint32_t bmp_offset;
  };
  
  struct bmpfile_info_header
  {
    uint32_t header_sz;
    int32_t width;
    int32_t height;
    uint16_t nplanes;
    uint16_t bitspp;
    uint32_t compress_type;
    uint32_t bmp_bytesz;
    int32_t hres;
    int32_t vres;
    uint32_t ncolors;
    uint32_t nimpcolors;
  };
  
  struct bmpfile
  {
    bmpfile_magic magic;
    bmpfile_header header;
    bmpfile_info_header info_header;
    char *data;
  };
}

#endif
