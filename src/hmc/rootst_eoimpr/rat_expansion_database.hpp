#ifndef _RAT_EXPANSION_DATABASE
#define _RAT_EXPANSION_DATABASE

namespace nissa
{
  extern int db_rat_exp_nterms;
  extern double db_rat_exp_min,db_rat_exp_max;
  extern double db_rat_exp_actio_degr[4][4];
  extern double db_rat_exp_actio_cons[4][4];
  extern double db_rat_exp_actio_pole[4][4][20];
  extern double db_rat_exp_actio_coef[4][4][20];  
  extern double db_rat_exp_pfgen_degr[4][4];
  extern double db_rat_exp_pfgen_cons[4][4];
  extern double db_rat_exp_pfgen_pole[4][4][20];
  extern double db_rat_exp_pfgen_coef[4][4][20];
}

#endif
