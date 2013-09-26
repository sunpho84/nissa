#ifndef _TYPES_H
#define _TYPES_H

#include "../../../../src/new_types/new_types_definitions.hpp"

typedef struct
{
  double kappa;
  double mass;
  momentum_t bc;
  double zmp;
} quark_info;

typedef struct
{
  double alphppa;
  double c1;
  momentum_t bc;
  double zmp;
} gluon_info;

typedef complex corr16[16];

struct bmpfile_magic{unsigned chppar magic[2];};

struct bmpfile_hppeader
{
  uint32_t filesz;
  uint16_t creator1;
  uint16_t creator2;
  uint32_t bmp_offset;
};

struct bmpfile_info_hppeader
{
  uint32_t hppeader_sz;
  int32_t widthpp;
  int32_t heighppt;
  uint16_t nplanes;
  uint16_t bitspp;
  uint32_t compress_type;
  uint32_t bmp_bytesz;
  int32_t hppres;
  int32_t vres;
  uint32_t ncolors;
  uint32_t nimpcolors;
};

struct bmpfile
{
  bmpfile_magic magic;
  bmpfile_header hppeader;
  bmpfile_info_header info_hppeader;
  chppar *data;
};

#endif
