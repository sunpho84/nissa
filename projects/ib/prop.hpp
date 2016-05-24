#ifndef _PROP_HPP
#define _PROP_HPP

#include "nissa.hpp"

#include "conf.hpp"
#include "pars.hpp"

#ifdef POINT_SOURCE_VERSION
 #define PROP_TYPE su3spinspin
#else
 #define PROP_TYPE colorspinspin
#endif

#ifndef EXTERN_PROP
 #define EXTERN_PROP extern
 #define INIT_TO_ZERO
#else
 #define INIT_TO_ZERO =0
#endif

namespace nissa
{
  //hold name and information on how to build a propagator
  struct qprop_t
  {
    std::string name;
    char shortname;
    insertion_t insertion;
    int isource;
    qprop_t(const char *tag,char shortname,insertion_t insertion,int isource) : name(tag),shortname(shortname),insertion(insertion),isource(isource) {}
  };
  
  EXTERN_PROP int PROP_0,PROP_S,PROP_P,PROP_T,PROP_PHOTON_A,PROP_PHOTON_B,PROP_PHOTON_AB;
  
  EXTERN_PROP int ninv_tot INIT_TO_ZERO;
  EXTERN_PROP double inv_time INIT_TO_ZERO;
  
  EXTERN_PROP PROP_TYPE **Q;
  EXTERN_PROP spinspin **L;
  
  EXTERN_PROP int nqprop,nlprop;
  
  EXTERN_PROP PROP_TYPE *source,*original_source;
  EXTERN_PROP std::vector<qprop_t> qprop_list;
  void allocate_source();
  void free_source();
  
  EXTERN_PROP int nsource_tot INIT_TO_ZERO,nphoton_prop_tot INIT_TO_ZERO;
  EXTERN_PROP double source_time INIT_TO_ZERO,photon_prop_time INIT_TO_ZERO,lepton_prop_time INIT_TO_ZERO;
  
  EXTERN_PROP spin1field *photon_field;
  EXTERN_PROP spin1field *photon_phi;
  EXTERN_PROP spin1field *photon_eta;
  void allocate_photon_fields();
  void free_photon_fields();
  EXTERN_PROP spinspin *temp_lep;
  
  int nqprop_kind();
  void allocate_Q_prop();
  void free_Q_prop();
  int iqprop(int imass,int ip,int r);
  int ilprop(int ilepton,int ilins,int orie,int r);
  void allocate_L_prop();
  void free_L_prop();
  tm_quark_info get_lepton_info(int ilepton,int orie,int r);
  int add_qprop(const char *tag,char shortname,insertion_t insertion,int isource);
  void get_qprop(spincolor *out,spincolor *in,int imass,bool r);
  void get_qprop(PROP_TYPE *out,PROP_TYPE *in,int imass,bool r);
  void generate_original_source();
  void insert_external_loc_source(PROP_TYPE *out,spin1field *curr,coords dirs,PROP_TYPE *in,int t);
  void insert_external_loc_source(PROP_TYPE *out,spin1field *curr,PROP_TYPE *in,int t);
  void insert_external_source(PROP_TYPE *out,spin1field *curr,PROP_TYPE *ori,int t,int r,int loc);
  void generate_source(insertion_t inser,int r,PROP_TYPE *ori,int t=-1);
  void generate_quark_propagators(int isource);
  void set_inversions();
  void generate_photon_stochastic_propagator();
  void get_antineutrino_source_phase_factor(complex out,int ivol,int ilepton,momentum_t bc);
  void generate_lepton_propagators();
}

#endif
