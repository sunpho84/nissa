#ifndef _TMLQCD_TYPES_H
#define _TMLQCD_TYPES_H

typedef struct
{
  _Complex double c0,c1,c2;
} su3_vector;

typedef struct
{
  su3_vector s0,s1,s2,s3;
} spinor;

typedef struct 
{
  _Complex double c00, c01, c02, c10, c11, c12, c20, c21, c22;
} tmlQCD_su3;

#endif
