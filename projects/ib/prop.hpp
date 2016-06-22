#ifndef _PROP_HPP
#define _PROP_HPP

#include "nissa.hpp"

#include "conf.hpp"
#include "pars.hpp"

#ifndef EXTERN_PROP
 #define EXTERN_PROP extern
#define INIT_TO(A)
#else
 #define INIT_TO(A) =A
#endif

namespace nissa
{
  inline int so_sp_col_ind(int sp,int col){return col+nso_col*sp;}
  
  //hold name and information on how to build a propagator
  struct qprop_t
  {
    bool is_source;
    
    double kappa;
    double mass;
    int r;
    double theta;
    insertion_t insertion;
    std::string source_name;
    int tins;
    double residue;
    
    rnd_t noise_type;
    double ori_source_norm2;
    
    //spincolor data
    std::vector<spincolor*> sp;
    spincolor* &operator[](int i){return sp[i];}
    void alloc_spincolor()
    {
      sp.resize(nso_spi*nso_col);
      for(int i=0;i<nso_spi*nso_col;i++)
	sp[i]=nissa_malloc("sp",loc_vol,spincolor);
  }
    
    //initialize as a propagator
    void init_as_propagator(insertion_t _insertion,std::string _source_name,int _tins,double _residue,double _kappa,double _mass,int _r,double _theta)
    {
      is_source=false;
      
      kappa=_kappa;
      mass=_mass;
      r=_r;
      theta=_theta;
      insertion=_insertion;
      source_name=_source_name;
      tins=_tins;
      residue=_residue;
      
      alloc_spincolor();
    }
    
    //initialize as a source
    void init_as_source(rnd_t _noise_type)
    {
      is_source=true;
      
      noise_type=_noise_type;
      alloc_spincolor();
    }
    
    qprop_t(insertion_t insertion,std::string source_name,int tins,double residue,double kappa,double mass,int r,double theta)
    {init_as_propagator(insertion,source_name,tins,residue,kappa,mass,r,theta);}
    qprop_t(rnd_t noise_type) {init_as_source(noise_type);}
    qprop_t() {is_source=0;}
    ~qprop_t() {for(size_t i=0;i<sp.size();i++) nissa_free(sp[i]);}
  };
  
  const int ALL_TIMES=-1;
  
  EXTERN_PROP int ninv_tot INIT_TO(0);
  EXTERN_PROP double inv_time INIT_TO(0);
  
  EXTERN_PROP int nstore_prop INIT_TO(0);
  EXTERN_PROP double store_prop_time INIT_TO(0);
  
  EXTERN_PROP int nread_prop INIT_TO(0);
  EXTERN_PROP double read_prop_time INIT_TO(0);
  
  EXTERN_PROP std::map<std::string,qprop_t> Q;
  EXTERN_PROP std::vector<std::string> qprop_name_list;
  EXTERN_PROP spinspin **L;
  
  EXTERN_PROP int nlprop;
  
  EXTERN_PROP std::vector<std::string> ori_source_name_list;
  EXTERN_PROP spincolor *loop_source;
  inline void allocate_loop_source(){loop_source=nissa_malloc("loop_source",loc_vol,spincolor);}
  inline void free_loop_source(){nissa_free(loop_source);}
  
  EXTERN_PROP int nsource_tot INIT_TO(0),nphoton_prop_tot INIT_TO(0);
  EXTERN_PROP double source_time INIT_TO(0),photon_prop_time INIT_TO(0),lepton_prop_time INIT_TO(0);
  
  EXTERN_PROP spin1field *photon_field;
  EXTERN_PROP spin1field *photon_phi;
  EXTERN_PROP spin1field *photon_eta;
  void allocate_photon_fields();
  void free_photon_fields();
  EXTERN_PROP spinspin *temp_lep;
  
  inline int ilprop(int ilepton,int ilins,int orie,int r)
  {return r+nr_lep*(ilins+nlins*(orie+norie*ilepton));}
  
  void allocate_L_prop();
  void free_L_prop();
  tm_quark_info get_lepton_info(int ilepton,int orie,int r);
  
  void get_qprop(spincolor *out,spincolor *in,double kappa,double mass,int r,double residue,double theta);
  void generate_original_source(qprop_t *sou);
  inline void generate_original_sources()
  {
    for(size_t i=0;i<ori_source_name_list.size();i++)
      generate_original_source(&Q[ori_source_name_list[i]]);
  }
  void insert_external_loc_source(spincolor *out,spin1field *curr,coords dirs,spincolor *in,int t);
  void insert_external_loc_source(spincolor *out,spin1field *curr,spincolor *in,int t);
  void insert_external_source(spincolor *out,spin1field *curr,spincolor *ori,int t,int r,int loc);
  void generate_source(insertion_t inser,int r,spincolor *ori,int t);
  void generate_quark_propagators(int isource);
  void generate_photon_stochastic_propagator();
  void get_antineutrino_source_phase_factor(complex out,int ivol,int ilepton,momentum_t bc);
  void generate_lepton_propagators();
  
  inline void start_hit(int ihit)
  {
    master_printf("\n=== Hit %d/%d ====\n",ihit+1,nhits);
    random_shift_gauge_conf(conf,old_theta,put_theta);
  }
  
  inline void generate_propagators(int ihit)
  {
    generate_photon_stochastic_propagator();
    generate_original_sources();
    generate_lepton_propagators();
    generate_quark_propagators(ihit);
  }
}

#undef INIT_TO

#endif
