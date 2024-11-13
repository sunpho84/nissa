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
  
  void clearCorrelations();
  
  void compute_prop_scalprod(complex& res,
			     const std::string& pr_dag,
			     const std::string& pr);
  
  ///////////////////////////////////////// meson contractions ///////////////////////////////////////////////////////
  
  EXTERN_CONTR std::string mes2pts_prefix INIT_TO("mes_contr");
  
  struct mes_contr_map_t
  {
    std::string name;
    std::string a,b;
    mes_contr_map_t(std::string name,std::string a,std::string b) : name(name),a(a),b(b) {}
  };
  EXTERN_CONTR std::vector<mes_contr_map_t> mes2pts_contr_map;
  EXTERN_CONTR int nmes2pts_contr_made INIT_TO(0);
  EXTERN_CONTR double mes2pts_contr_time INIT_TO(0);
  CUDA_MANAGED EXTERN_CONTR LxField<complex> *loc_contr;
  
  CUDA_MANAGED EXTERN_CONTR complex *mes2pts_contr INIT_TO(NULL);
  EXTERN_CONTR std::vector<idirac_pair_t> mes_gamma_list;
  void allocate_mes2pts_contr();
  void compute_mes2pt_contr(int icombo);
  void print_mes2pts_contr(int n=nhits,int force_append=false,int skip_inner_header=false,const std::string &alternative_header_template="");
  void free_mes2pts_contr();
  
  inline int ind_mes2pts_contr(int iquark_combo,int ihadr_contr,int t)
  {
    return
      (t+glbSize[0]*
       (ihadr_contr+mes_gamma_list.size()*
	iquark_combo));
  }
  EXTERN_CONTR int mes2pts_contr_size;
  
  ///////////////////////////////////////// handcuffs contractions ///////////////////////////////////////////////////////
  
  void vector_current_mel(LxField<spin1field>& si,
			  const dirac_matr& ext_g,
			  const int& r,
			  const char* id_Qbw,
			  const char* id_Qfw,
			  const bool& revert);
  
  void conserved_vector_current_mel(const LxField<quad_su3>& conf,
				    LxField<spin1field>& si,
				    const dirac_matr& ext_g,
				    const int& r,
				    const char* name_bw,
				    const char* name_fw,
				    const bool& revert);
  
  void local_or_conserved_vector_current_mel(LxField<spin1field>& si,
					     const dirac_matr &g,
					     const std::string &prop_name_bw,
					     const std::string &prop_name_fw,
					     const bool& revert);
  
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
    
    std::string left;
    
    std::string right;
    
    handcuffs_map_t(const std::string& name,
		    const std::string& left,
		    const std::string& right) :
      name(name),
      left(left),
      right(right)
    {
    }
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
  CUDA_MANAGED EXTERN_CONTR dirac_matr Cg5;
  void set_Cg5();
  
  EXTERN_CONTR int compute_octet INIT_TO(0);
  EXTERN_CONTR int compute_decuplet INIT_TO(0);
  EXTERN_CONTR std::vector<bar_triplet_t> bar2pts_contr_map;
  EXTERN_CONTR int nbar2pts_contr_made INIT_TO(0);
  EXTERN_CONTR double bar2pts_contr_time INIT_TO(0);
  CUDA_MANAGED EXTERN_CONTR complex *bar2pts_contr INIT_TO(NULL);
  void set_bar2pts_contr_ins_map();
  void allocate_bar2pts_contr();
  void compute_bar2pts_contr();
  void print_bar2pts_contr();
  void free_bar2pts_contr();
  
  inline int ind_bar2pts_contr(int icombo,int dir_exc,int t)
  {
    return
      (t+glbSize[0]*
       (dir_exc+2*
	icombo));
  }
  EXTERN_CONTR int bar2pts_contr_size;
  
#define BAR_ALT_LIMITED_PROJ

#ifdef BAR_ALT_LIMITED_PROJ
 #define NBAR_ALT_PROJ 3
 #define NBAR_ALT_SINGLE_PROJ 10
#else
 #define NBAR_ALT_PROJ 6
 #define NBAR_ALT_SINGLE_PROJ 17
#endif

  CUDA_MANAGED EXTERN_CONTR complex *bar2pts_alt_contr INIT_TO(NULL);
  EXTERN_CONTR int nbar2pts_alt_contr_made INIT_TO(0);
  EXTERN_CONTR double bar2pts_alt_contr_time INIT_TO(0);
  void allocate_bar2pts_alt_contr();
  void compute_bar2pts_alt_contr();
  void print_bar2pts_alt_contr();
  void free_bar2pts_alt_contr();
  inline int ind_bar2pts_alt_contr(int icombo,int iWick,int iProj,int t)
  {
    return
      (t+glbSize[0]*
       (iProj+NBAR_ALT_PROJ*
	(iWick+2*
	 icombo)));
  }
  EXTERN_CONTR int bar2pts_alt_contr_size;
  
  EXTERN_CONTR int nsmear_oper INIT_TO(0);
  EXTERN_CONTR double smear_oper_time INIT_TO(0);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //compute all contractions
  inline void compute_contractions()
  {
    compute_handcuffs_contr();
    //compute_meslep_contr();
    if(compute_octet) crash("dependencies are broken");//compute_bar2pts_contr();
    if(compute_decuplet) crash("dependencies are broken");//compute_bar2pts_alt_contr();
  }
  
  //print out all contractions
  inline void print_contractions()
  {
    print_mes2pts_contr();
    print_handcuffs_contr();
    //print_meslep_contr();
    if(compute_octet) print_bar2pts_contr();
    if(compute_decuplet) print_bar2pts_alt_contr();
  }
}

#undef INIT_TO

#endif
