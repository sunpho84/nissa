#ifndef _CONTR_HPP
#define _CONTR_HPP

#include <vector>

#include "geometry/geometry_lx.hpp"
#include "new_types/dirac.hpp"

#include "pars.hpp"

#ifndef EXTERN_CONTR
 #define EXTERN_CONTR extern
 #define INIT_TO(VAR)
#else
 #define INIT_TO(VAR) =VAR
#endif

namespace nissa
{
  EXTERN_CONTR double contr_print_time INIT_TO(0);
  
  ///////////////////////////////////////// meson contractions ///////////////////////////////////////////////////////
  
  struct mes_doublet_t
  {
    int a,b;
    mes_doublet_t(int a,int b) : a(a),b(b) {}
  };
  EXTERN_CONTR std::vector<mes_doublet_t> mes2pts_contr_ins_map;
  EXTERN_CONTR std::vector<mes_doublet_t> mes2pts_contr_quark_map;
  EXTERN_CONTR int nmes2pts_contr INIT_TO(0);
  EXTERN_CONTR double mes2pts_contr_time INIT_TO(0);
  EXTERN_CONTR complex *mes2pts_contr INIT_TO(NULL);
  EXTERN_CONTR std::vector<idirac_pair_t> mes_gamma_list;
  void set_mes2pts_contr_ins_map();
  void allocate_mes2pts_contr();
  void compute_mes2pts_contr();
  void print_mes2pts_contr();
  void free_mes2pts_contr();
  
  inline int ind_mes2pts_contr(int ins,int iquark_combo,int ihadr_contr,int t)
  {
    return
      (t+glb_size[0]*
       (ihadr_contr+mes_gamma_list.size()*
	(iquark_combo+mes2pts_contr_quark_map.size()*
	 ins)));
  }
  EXTERN_CONTR int mes2pts_contr_size;
  
  ////////////////////////////////////  mesoleptonic contraction /////////////////////////////////////////////////////
  
  EXTERN_CONTR int nmeslep_contr INIT_TO(0);
  EXTERN_CONTR double meslep_contr_time INIT_TO(0);

  //list the 8 matrices to insert for the weak current
  const int nmeslep_weak_ins=17;
  const int nindep_meslep_weak=9;
  EXTERN_CONTR int nmeslep_corr;
  EXTERN_CONTR spinspin *meslep_hadr_part;
  EXTERN_CONTR complex *meslep_contr;
  
  //parameters of the leptons
  EXTERN_CONTR int *lep_contr_iq1;
  EXTERN_CONTR int *lep_contr_iq2;
  
  //const int nmeslep_proj=4,meslep_projs[nmeslep_proj]={9,4,5,0};
  const int nmeslep_proj=1,meslep_projs[nmeslep_proj]={4};
  const int list_weak_insq[nmeslep_weak_ins]=     {1,2,3,4, 6,7,8,9,  1,2,3,4, 6,7,8,9, 5};
  const int list_weak_insl[nmeslep_weak_ins]=     {1,2,3,4, 6,7,8,9,  6,7,8,9, 1,2,3,4, 5};
  const int list_weak_ind_contr[nmeslep_weak_ins]={0,0,0,1, 2,2,2,3,  4,4,4,5, 6,6,6,7, 8};
  const char list_weak_ind_nameq[nindep_meslep_weak][3]={"VK","V0","AK","A0","VK","V0","AK","A0","P5"};
  const char list_weak_ind_namel[nindep_meslep_weak][3]={"VK","V0","AK","A0","AK","A0","VK","V0","V0"};
  
  void compute_meslep_contr();
  void print_meslep_contr();
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //barion contractions
  struct bar_triplet_t
  {
    int a,b,c;
    bar_triplet_t(int a,int b,int c) : a(a),b(b),c(c) {}
  };
  
  //Cg5
  EXTERN_CONTR dirac_matr Cg5;
  void set_Cg5();
  
  EXTERN_CONTR std::vector<bar_triplet_t> prop_bar_contr_map;
  EXTERN_CONTR int nbar_contr INIT_TO(0);
  EXTERN_CONTR double bar_contr_time INIT_TO(0);
  EXTERN_CONTR complex *bar_contr INIT_TO(NULL);
  void set_bar_prop_contr_list();
  void allocate_bar_contr();
  void compute_bar_contr();
  void print_bar_contr();
  void free_bar_contr();
  
  inline int ind_bar_contr(int icombo,int ism_sink,int iqa,int iqb,int iqc,int dir_exc,int t)
  {return
      (t+glb_size[0]*
       (dir_exc+2*
	(iqc+nquarks*
	 (iqb+nquarks*
	  (iqa+nquarks*
	   (ism_sink+nsm_sink*icombo))))));
  }
  EXTERN_CONTR int bar_contr_size;
  
  EXTERN_CONTR int nsmear_oper INIT_TO(0);
  EXTERN_CONTR double smear_oper_time INIT_TO(0);
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //compute all contractions
  inline void compute_contractions()
  {
    if(compute_mes2pts_flag) compute_mes2pts_contr();
    if(compute_meslep_flag) compute_meslep_contr();
    if(compute_bar_flag) compute_bar_contr();
  }
  
  //print out all contractions
  inline void print_contractions()
  {
    if(compute_mes2pts_flag) print_mes2pts_contr();
    if(compute_meslep_flag) print_meslep_contr();
    if(compute_bar_flag) print_bar_contr();
  }
}

#undef INIT_TO

#endif
