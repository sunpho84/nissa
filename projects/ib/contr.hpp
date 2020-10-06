#ifndef _CONTR_HPP
#define _CONTR_HPP

#include <vector>

#include "geometry/geometry_lx.hpp"
#include "new_types/dirac.hpp"

#include "pars.hpp"

#include "meslep.hpp"


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
  
  struct mes_contr_map_t
  {
    std::string name;
    std::string a,b;
    mes_contr_map_t(std::string name,std::string a,std::string b) : name(name),a(a),b(b) {}
  };
  EXTERN_CONTR std::vector<mes_contr_map_t> mes2pts_contr_map;
  EXTERN_CONTR int nmes2pts_contr_made INIT_TO(0);
  EXTERN_CONTR double mes2pts_contr_time INIT_TO(0);
  EXTERN_CONTR complex *mes2pts_contr INIT_TO(NULL);
  EXTERN_CONTR std::vector<idirac_pair_t> mes_gamma_list;
  void allocate_mes2pts_contr();
  void compute_mes2pts_contr();
  void print_mes2pts_contr();
  void free_mes2pts_contr();
  
  inline int ind_mes2pts_contr(int iquark_combo,int ihadr_contr,int t)
  {
    return
      (t+glb_size[0]*
       (ihadr_contr+mes_gamma_list.size()*
	iquark_combo));
  }
  EXTERN_CONTR int mes2pts_contr_size;
  
  ///////////////////////////////////////// handcuffs contractions ///////////////////////////////////////////////////////
  
  struct handcuffs_side_map_t
  {
    std::string name;
    int igamma;
    std::string bw,fw;
    int store;
    handcuffs_side_map_t(std::string name,int igamma,std::string bw,std::string fw,int store) : name(name),igamma(igamma),bw(bw),fw(fw),store(store) {}
  };
  EXTERN_CONTR std::vector<handcuffs_side_map_t> handcuffs_side_map;
  
  struct handcuffs_map_t
  {
    std::string name;
    std::string left,right;
    handcuffs_map_t(std::string name,std::string left,std::string right) : name(name),left(left),right(right) {}
  };
  EXTERN_CONTR std::vector<handcuffs_map_t> handcuffs_map;
  EXTERN_CONTR int nhandcuffs_contr_made INIT_TO(0);
  EXTERN_CONTR double handcuffs_contr_time INIT_TO(0);
  EXTERN_CONTR complex *handcuffs_contr INIT_TO(NULL);
  void allocate_handcuffs_contr();
  void compute_handcuffs_contr();
  void print_handcuffs_contr();
  void free_handcuffs_contr();
  
  inline int ind_handcuffs_contr(int ihand)
  {return ihand;}
  EXTERN_CONTR int handcuffs_contr_size;
  
  //////////////////////////////////////////////// barion contractions //////////////////////////////////////////////////
  
  struct bar_triplet_t
  {
    std::string name;
    std::string a,b,c;
    bar_triplet_t(std::string name,std::string a,std::string b,std::string c) : name(name),a(a),b(b),c(c) {}
  };
  
  //Cg5
  EXTERN_CONTR dirac_matr Cg5;
  void set_Cg5();
  
  EXTERN_CONTR std::vector<bar_triplet_t> bar2pts_contr_map;
  EXTERN_CONTR int nbar2pts_contr_made INIT_TO(0);
  EXTERN_CONTR double bar2pts_contr_time INIT_TO(0);
  EXTERN_CONTR complex *bar2pts_contr INIT_TO(NULL);
  void set_bar2pts_contr_ins_map();
  void allocate_bar2pts_contr();
  void compute_bar2pts_contr();
  void print_bar2pts_contr();
  void free_bar2pts_contr();
  
  inline int ind_bar2pts_contr(int icombo,int dir_exc,int t)
  {
    return
      (t+glb_size[0]*
       (dir_exc+2*
	(icombo+bar2pts_contr_map.size())));
  }
  EXTERN_CONTR int bar2pts_contr_size;
  
  EXTERN_CONTR int nsmear_oper INIT_TO(0);
  EXTERN_CONTR double smear_oper_time INIT_TO(0);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //compute all contractions
  inline void compute_contractions()
  {
    compute_mes2pts_contr();
    compute_handcuffs_contr();
    //compute_meslep_contr();
    compute_bar2pts_contr();
  }
  
  //print out all contractions
  inline void print_contractions()
  {
    print_mes2pts_contr();
    print_handcuffs_contr();
    //print_meslep_contr();
    print_bar2pts_contr();
  }
}

#undef INIT_TO

#endif
