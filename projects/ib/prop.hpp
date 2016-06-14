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
  //hold name and information on how to build a propagator
  struct qprop_t
  {
    std::string name;
    char shortname;
    insertion_t insertion;
    int isource;
    int tins;
    qprop_t(const char *tag,char shortname,insertion_t insertion,int isource,int tins) :
      name(tag),shortname(shortname),insertion(insertion),isource(isource),tins(tins) {}
  };
  
  const int ALL_TIMES=-1;
  EXTERN_PROP int ORI_SOURCE INIT_TO(-1);
  EXTERN_PROP int PROP_0,PROP_S,PROP_P,PROP_T,PROP_PHOTON_A,PROP_PHOTON_B,PROP_PHOTON_AB;
  
  EXTERN_PROP int ninv_tot INIT_TO(0);
  EXTERN_PROP double inv_time INIT_TO(0);
  
  EXTERN_PROP spincolor **Q;
  EXTERN_PROP spinspin **L;
  
  EXTERN_PROP int nqprop,nlprop;
  
  EXTERN_PROP spincolor *source;
  EXTERN_PROP std::vector<qprop_t> qprop_list;
  void allocate_source();
  void free_source();
  
  EXTERN_PROP int nsource_tot INIT_TO(0),nphoton_prop_tot INIT_TO(0);
  EXTERN_PROP double source_time INIT_TO(0),photon_prop_time INIT_TO(0),lepton_prop_time INIT_TO(0);
  
  EXTERN_PROP spin1field *photon_field;
  EXTERN_PROP spin1field *photon_phi;
  EXTERN_PROP spin1field *photon_eta;
  void allocate_photon_fields();
  void free_photon_fields();
  EXTERN_PROP spinspin *temp_lep;
  
  int nqprop_kind();
  void allocate_Q_prop();
  void free_Q_prop();
  
  inline int iqprop(int iquark,int ip,int is,int ic)
  {
    int c;
    if(ip==ORI_SOURCE) c=0;
    else c=1+iquark+nquarks*ip;
    
    return ic+nso_col*(is+nso_spi*c);
  }
  
  inline int ilprop(int ilepton,int ilins,int orie,int r)
  {return r+nr_lep*(ilins+nlins*(orie+norie*ilepton));}
  
  void allocate_L_prop();
  void free_L_prop();
  tm_quark_info get_lepton_info(int ilepton,int orie,int r);
  int add_qprop(const char *tag,char shortname,insertion_t insertion,int isource,int tins=ALL_TIMES);
  void get_qprop(spincolor *out,spincolor *in,int iquark);
  void generate_original_source();
  void insert_external_loc_source(spincolor *out,spin1field *curr,coords dirs,spincolor *in,int t);
  void insert_external_loc_source(spincolor *out,spin1field *curr,spincolor *in,int t);
  void insert_external_source(spincolor *out,spin1field *curr,spincolor *ori,int t,int r,int loc);
  void generate_source(insertion_t inser,int r,spincolor *ori,int t);
  void generate_quark_propagators(int isource);
  void set_inversions();
  void generate_photon_stochastic_propagator();
  void get_antineutrino_source_phase_factor(complex out,int ivol,int ilepton,momentum_t bc);
  void generate_lepton_propagators();
  
  inline void start_source(int isource)
  {
    master_printf("\n=== Source %d/%d ====\n",isource+1,nsources);
    random_shift_gauge_conf(conf,old_theta,put_theta);
  }
  
  inline void generate_propagators(int isource)
  {
    generate_photon_stochastic_propagator();
    generate_original_source();
    if(compute_meslep_flag) generate_lepton_propagators();
    generate_quark_propagators(isource);
  }
}

#endif
