#ifndef _MESLEP_HPP
#define _MESLEP_HPP

#include "pars.hpp"

#ifndef EXTERN_MESLEP
 #define EXTERN_MESLEP extern
 #define INIT_TO(VAR)
#else
 #define INIT_TO(VAR) =VAR
#endif

namespace nissa
{
  EXTERN_MESLEP int nquark_lep_combos;
  EXTERN_MESLEP double *lep_energy,*neu_energy;
  
  EXTERN_MESLEP int nmeslep_contr_made INIT_TO(0);
  EXTERN_MESLEP double meslep_contr_time INIT_TO(0);
  
  //sign of the lepton momentum
  const int sign_orie[2]={-1,+1};
  
  EXTERN_MESLEP tm_quark_info *leps;
  
  //list the 8 matrices to insert for the weak current
  const int nmeslep_weak_ins=17;
  const int nindep_meslep_weak=9;
  EXTERN_MESLEP int nmeslep_corr;
  EXTERN_MESLEP spinspin *meslep_hadr_part;
  EXTERN_MESLEP complex *meslep_contr;
  
  //parameters of the leptons
  EXTERN_MESLEP int *lep_contr_iq1;
  EXTERN_MESLEP int *lep_contr_iq2;
  
  //const int nmeslep_proj=4,meslep_projs[nmeslep_proj]={9,4,5,0};
  const int nmeslep_proj=1,meslep_projs[nmeslep_proj]={4};
  const int list_weak_insq[nmeslep_weak_ins]=     {1,2,3,4, 6,7,8,9,  1,2,3,4, 6,7,8,9, 5};
  const int list_weak_insl[nmeslep_weak_ins]=     {1,2,3,4, 6,7,8,9,  6,7,8,9, 1,2,3,4, 5};
  const int list_weak_ind_contr[nmeslep_weak_ins]={0,0,0,1, 2,2,2,3,  4,4,4,5, 6,6,6,7, 8};
  const char list_weak_ind_nameq[nindep_meslep_weak][3]={"VK","V0","AK","A0","VK","V0","AK","A0","P5"};
  const char list_weak_ind_namel[nindep_meslep_weak][3]={"VK","V0","AK","A0","AK","A0","VK","V0","V0"};
  
  void allocate_L_prop();
  void free_L_prop();
  tm_quark_info get_lepton_info(int ilepton,int orie,int r);
  
  void read_meslep_contr_pars();
  void compute_meslep_contr();
  void print_meslep_contr();
}

#undef INIT_TO

#endif
