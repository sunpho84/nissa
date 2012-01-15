#pragma once

typedef struct
{
  char name[20];
  double con;
  int npoles;
  double *poles;
  double *weights;
} rat_approx;

typedef struct
{
  int deg;
  double mass;
  double imm_pot;
} quark_content;

