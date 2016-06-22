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
  
  
  struct tm_quark_info
  {
    double kappa;
    double mass;
    momentum_t bc;
    double zmp;
    int r;
    tm_quark_info(double kappa,double mass,int r,double theta) :
      kappa(kappa),mass(mass),r(r) {bc[0]=1;for(int mu=1;mu<NDIM;mu++) bc[mu]=theta;}
    tm_quark_info() {}
  };
  
  struct gauge_info
  {
    zero_mode_sub_t zms;
    double alpha;
    double c1;
    momentum_t bc;
    gauge_info()
    {
      zms=UNNO_ALEMANNA;
      alpha=LANDAU_ALPHA;
      c1=WILSON_C1;
      for(int mu=0;mu<NDIM;mu++) bc[mu]=0;
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
