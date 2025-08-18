#ifndef _CONTR_HPP
#define _CONTR_HPP

#include <complex>
#include <vector>

#include "geometry/geometry_lx.hpp"
#include "ib/prop.hpp"
#include "new_types/dirac.hpp"

#include "pars.hpp"

#include "meslep.hpp"


#ifndef EXTERN_CONTR
# define EXTERN_CONTR extern
# define INIT_TO(VAR)
#else
# define INIT_TO(VAR) =VAR
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
  
  struct mes_contr_t
  {
    std::string name;
    
    std::string a;
    
    std::string b;
    
    std::vector<idirac_pair_t> gammaList;
    
    complex* contr;
    
    size_t contrSize;
    
    /// Allocate mesonic contractions
    void alloc()
    {
      contrSize=glbSize[0]*gammaList.size();
      contr=new complex[contrSize];
    }
    
    /// Allocate mesonic contractions
    void unalloc()
    {
      delete[] contr;
    }
    
    complex& operator()(const size_t& iGamma,
			const size_t& t)
    {
      return contr[t+glbSize[0]*iGamma];
    }
    
    /// Reset the contractions
    void reset()
    {
      for(size_t i=0;i<contrSize;i++)
	complex_put_to_zero(contr[i]);
    }
    
    mes_contr_t(const std::string& name,
		const std::string& a,
		const std::string& b) :
      name(name),
      a(a),
      b(b)
    {
    }
  };
  
  EXTERN_CONTR std::vector<mes_contr_t> mes2ptsContr;
  
  EXTERN_CONTR int nmes2pts_contr_made INIT_TO(0);
  EXTERN_CONTR double mes2pts_contr_time INIT_TO(0);
  
  EXTERN_CONTR int nmes2pts_move_to_make_readable_made INIT_TO(0);
  EXTERN_CONTR double mes2pts_move_to_make_readable_time INIT_TO(0);
  
  CUDA_MANAGED EXTERN_CONTR LxField<complex> *loc_contr;
  
  using ContrProp=LxField<spincolor,defaultSpaceTimeLayout,defaultMemorySpace>;
  EXTERN_CONTR std::map<std::string,std::vector<ContrProp*>> mes2ptsPropsLib;
  inline void removeMes2PtsProp(const std::string& n)
  {
    MASTER_PRINTF("Removing %s from the mes2ptsPropsLib\n",n.c_str());
    if constexpr(defaultMemorySpace!=MemorySpace::CPU)
      for(auto& vi : mes2ptsPropsLib[n])
	delete vi;
    mes2ptsPropsLib.erase(n);
  }
  
  void compute_mes2pt_contr(const size_t& icombo);
  
  void print_mes2pts_contr(const int iHit,int n=nHits,int force_append=false,int skip_inner_header=false,const std::string &alternative_header_template="");
  void free_mes2pts_contr();
  
  inline int ind_mes2pts_contr(const int& ihadr_contr,
			       const int& t)
  {
    return
      t+glbSize[0]*
      ihadr_contr;
  }
  
  ///////////////////////////////////////// handcuffs contractions ///////////////////////////////////////////////////////
  
  EXTERN_CONTR int computeHitSummedHandcuffs;
  
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
  
  /// Contains handcuffs
  struct HandcuffsSide
  {
    /// Gamma matrix
    int igamma;
    
    /// Name of the backward propagator
    std::string bw;
    
    /// Name of the forward propagator
    std::string fw;
    
    /// Decide whether to store or not
    int store;
    
    /// Handcuff data
    LxField<spin1field> data;
    
    /// Sum of the hit-by-hit handcuffs
    LxField<spin1field>* sum;
    
    /// Constructor
    HandcuffsSide(const int& igamma,
		  const std::string& bw,
		  const std::string& fw,
		  const int& store) :
      igamma(igamma),
      bw(bw),
      fw(fw),
      store(store),
      data("data"),
      sum(nullptr)
    {
      if(computeHitSummedHandcuffs)
	sum=new LxField<spin1field>("sum");
    }
    
    /// Move constructor
    HandcuffsSide(HandcuffsSide&&)=default;
    
    /// Destructor
    ~HandcuffsSide()
    {
      if(sum)
	delete sum;
    }
  };
  
  EXTERN_CONTR std::map<std::string,HandcuffsSide> handcuffsSides;
  
  /// Define a handcuff
  struct Handcuff
  {
    /// Name of the handcuff
    std::string name;
    
    /// Left side
    std::string left;
    
    /// Right side
    std::string right;
    
    /// Data of the handcuffs
    std::vector<std::complex<double>> data;
    
    /// Sum of the hit-by-hit handcuff
    std::complex<double> sum;
    
    /// Constructor
    Handcuff(const std::string& name,
	     const std::string& left,
	     const std::string& right) :
      name(name),
      left(left),
      right(right)
    {
    }
  };
  
  EXTERN_CONTR std::vector<Handcuff> handcuffs;
  
  EXTERN_CONTR double handcuffsContrTime INIT_TO(0);
  
  EXTERN_CONTR int nhandcuffsContrMade INIT_TO(0);
  
  void computeHandcuffSide(HandcuffsSide& h);
  
  void computeHandcuffsContr(const bool& useSum);
  
  void printHandcuffsContr(const int iHit);
  
  void freeHandcuffsContr();
  
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
  void print_bar2pts_contr(const int iHit);
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
  void print_bar2pts_alt_contr(const int iHit);
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
    //compute_meslep_contr();
    if(compute_octet) CRASH("dependencies are broken");//compute_bar2pts_contr();
    if(compute_decuplet) CRASH("dependencies are broken");//compute_bar2pts_alt_contr();
  }
  
  //print out all contractions
  inline void print_contractions(const int iHit=0)
  {
    print_mes2pts_contr(iHit);
    printHandcuffsContr(iHit);
    //print_meslep_contr();
    if(compute_octet) print_bar2pts_contr(iHit);
    if(compute_decuplet) print_bar2pts_alt_contr(iHit);
  }
}

#undef INIT_TO

#endif
