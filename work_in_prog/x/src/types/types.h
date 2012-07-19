#ifndef _TYPES_H
#define _TYPES_H

#include "../../../../src/new_types/new_types_definitions.h"

typedef struct
{
  double kappa;
  double mass;
  momentum_t bc;
  double zmp;
} quark_info;

typedef struct
{
  double alpha;
  double c1;
  momentum_t bc;
  double zmp;
} gluon_info;

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

#endif
