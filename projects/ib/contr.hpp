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
  //meson contractions
  struct mes_doublet_t
  {
    int a,b;
    mes_doublet_t(int a,int b) : a(a),b(b) {}
  };
  EXTERN_CONTR std::vector<mes_doublet_t> prop_mes_contr_map;
  EXTERN_CONTR int nmes_contract INIT_TO(0);
  EXTERN_CONTR double mes_contract_time INIT_TO(0);
  void set_mes_contract_list();
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //Cg5
  EXTERN_CONTR dirac_matr Cg5;
  void set_Cg5();
  
  //barion contractions
  struct bar_triplet_t
  {
    int a,b,c;
    bar_triplet_t(int a,int b,int c) : a(a),b(b),c(c) {}
  };
  EXTERN_CONTR std::vector<bar_triplet_t> prop_bar_contr_map;
  EXTERN_CONTR int nbar_contract INIT_TO(0);
  EXTERN_CONTR double bar_contract_time INIT_TO(0);
  EXTERN_CONTR complex *bar_contr INIT_TO(NULL);
  void set_bar_contract_list();
  void compute_bar_contractions();
  void print_bar_contractions();
  
  inline int ind_bar_contr(int icombo,int ism_sink,int ima,int ra,int imb,int rb,int imc,int rc,int dir_exc,int t)
  {return
      (t+glb_size[0]*
       (dir_exc+2*
	(rc+nr*
	 (imc+nqmass*
	  (rb+nr*
	   (imb+nqmass*
	    (ra+nr*
	     (ima+nqmass*
	      (ism_sink+nsm_sink*icombo)))))))));
  }
  EXTERN_CONTR int bar_contr_size;
  
  EXTERN_CONTR int nsmear_oper INIT_TO(0);
  EXTERN_CONTR double smear_oper_time INIT_TO(0);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //mesoleptonic contractions
  EXTERN_CONTR int nmeslep_contract INIT_TO(0);
  EXTERN_CONTR double meslep_contract_time INIT_TO(0);
}

#undef INIT_TO

#endif
