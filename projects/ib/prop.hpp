#ifndef _PROP_HPP
#define _PROP_HPP

#include "nissa.hpp"

#include "conf.hpp"
#include "meslep.hpp"
#include "pars.hpp"

#ifndef EXTERN_PROP
 #define EXTERN_PROP extern
#define INIT_TO(A)
#else
 #define INIT_TO(A) =A
#endif

namespace nissa
{
  //keep trace if generating photon is needed
  EXTERN_PROP int need_photon INIT_TO(0);

  CUDA_HOST_AND_DEVICE
  inline int so_sp_col_ind(const int& sp,const int& col)
  {
    return col+nso_col*sp;
  }
  
  typedef std::pair<std::string,std::pair<double,double>> source_term_t;
  
  //hold name and information on how to build a propagator
  struct qprop_t
  {
    bool is_source;
    
    double kappa;
    Momentum kappa_asymm;
    double mass;
    int r;
    double charge;
    Momentum theta;
    
    insertion_t insertion;
    std::vector<source_term_t> source_terms;
    GlbCoord tins;
    double residue;
    bool store;
    
    char ext_field_path[33];
    
    rnd_t noise_type;
    double ori_source_norm2;
    
    //spincolor data
    spincolor** sp;
    
    CUDA_HOST_AND_DEVICE
    spincolor* const &operator[](const int i) const
    {
      return sp[i];
    }
    
    CUDA_HOST_AND_DEVICE
    spincolor* &operator[](const int i)
    {
      return sp[i];
    }
    
    void alloc_spincolor()
    {
      sp=nissa_malloc("sp",nso_spi*nso_col,spincolor*);
      for(int i=0;i<nso_spi*nso_col;i++)
	sp[i]=nissa_malloc("sp",locVolWithBord.nastyConvert(),spincolor);
    }
    
    //initialize as a propagator
    void init_as_propagator(insertion_t _insertion,const std::vector<source_term_t>& _source_terms,const GlbCoord& _tins,double _residue,double _kappa,const Momentum& _kappa_asymm,double _mass,char *_ext_field_path,int _r,double _charge,const Momentum& _theta,bool _store)
    {
      is_source=false;
      
      kappa=_kappa;
      FOR_ALL_DIRS(mu)
	kappa_asymm(mu)=_kappa_asymm(mu);
      mass=_mass;
      r=_r;
      charge=_charge;
      FOR_ALL_DIRS(mu)
	theta(mu)=_theta(mu);
      insertion=_insertion;
      source_terms=_source_terms;
      tins=_tins;
      strncpy(ext_field_path,_ext_field_path,32);
      residue=_residue;
      store=_store;
      
      if(is_photon_ins(insertion)) need_photon=true;
      
      alloc_spincolor();
    }
    
    //initialize as a source
    void init_as_source(rnd_t _noise_type,const GlbCoord& _tins,int _r,bool _store)
    {
      is_source=true;
      
      noise_type=_noise_type;
      tins=_tins;
      r=_r;
      store=_store;
      alloc_spincolor();
    }
    
    qprop_t(insertion_t insertion,const std::vector<source_term_t>& source_terms,const GlbCoord& tins,double residue,double kappa,const Momentum& kappa_asymm,double mass,char *ext_field_path,int r,double charge,const Momentum& theta,bool store)
    {
      init_as_propagator(insertion,source_terms,tins,residue,kappa,kappa_asymm,mass,ext_field_path,r,charge,theta,store);
    }
    
    qprop_t(rnd_t noise_type,const GlbCoord& tins,int r,bool store)
    {
      init_as_source(noise_type,tins,r,store);
    }
    
    qprop_t()
    {
      is_source=0;
    }
    
    ~qprop_t()
    {
      for(int i=0;i<nso_spi*nso_col;i++)
	nissa_free(sp[i]);
      nissa_free(sp);
    }
  };
  
  const int ALL_TIMES=-1;
  
  EXTERN_PROP int ninv_tot INIT_TO(0);
  EXTERN_PROP double inv_time INIT_TO(0);
  
  EXTERN_PROP int nsme_tot INIT_TO(0);
  EXTERN_PROP double sme_time INIT_TO(0);
  
  EXTERN_PROP int nflw_tot INIT_TO(0);
  EXTERN_PROP double flw_time INIT_TO(0);
  
  EXTERN_PROP int nbflw_tot INIT_TO(0);
  EXTERN_PROP double bflw_time INIT_TO(0);
  
  EXTERN_PROP int nfft_tot INIT_TO(0);
  EXTERN_PROP double fft_time INIT_TO(0);
  void init_fft_filter_from_range(std::vector<std::pair<fft_mom_range_t,double>>& fft_mom_range_list);
  void init_fft_filterer_from_file(const char *fileout_suff,const char *filein_name);

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
  inline void allocate_loop_source()
  {
    loop_source=nissa_malloc("loop_source",locVolWithBord.nastyConvert(),spincolor);
  }
  
  inline void free_loop_source()
  {
    nissa_free(loop_source);
  }
  
  EXTERN_PROP int nsource_tot INIT_TO(0),nphoton_prop_tot INIT_TO(0);
  EXTERN_PROP double source_time INIT_TO(0),photon_prop_time INIT_TO(0),lepton_prop_time INIT_TO(0);
  
  EXTERN_PROP int load_photons INIT_TO(false);
  EXTERN_PROP int store_photons INIT_TO(false);
  CUDA_MANAGED EXTERN_PROP spin1field *photon_field;
  EXTERN_PROP spin1field *photon_phi;
  EXTERN_PROP spin1field *photon_eta;
  void allocate_photon_fields();
  void free_photon_fields();
  CUDA_MANAGED EXTERN_PROP spinspin *temp_lep;
  
  void get_qprop(spincolor *out,spincolor *in,double kappa,double mass,int r,double q,double residue,const Momentum& theta);
  void generate_original_source(qprop_t *sou);
  void generate_original_sources(int ihit,bool skip_io=false);
  void insert_external_loc_source(spincolor *out,spin1field *curr,spincolor *in,const GlbCoord& t,const Coords<bool>& dirs);
  void insert_external_source(spincolor *out,quad_su3 *conf,spin1field *curr,spincolor *ori,const GlbCoord& t,int r,const Coords<bool>& dirs,int loc);
  void generate_source(insertion_t inser,int r,double charge,double kappa,const Momentum& theta,spincolor *ori,const GlbCoord& t);
  void generate_quark_propagators(int isource);
  void generate_photon_stochastic_propagator(int ihit);
  //CUDA_HOST_AND_DEVICE void get_antineutrino_source_phase_factor(complex out,const int ivol,const int ilepton,const momentum_t bc);
  void generate_lepton_propagators();
  void propagators_fft(int ihit);
  
  void add_photon_field_to_conf(quad_su3 *conf,double charge);
  
  struct fft_filterer_t
  {
    int nfft_filtered;
    all_to_all_comm_t fft_filter_remap;
    std::string file_suffix;
    fft_filterer_t(int nfft_filtered,const all_to_all_scattering_list_t &sl,const std::string &file_suffix)
      : nfft_filtered(nfft_filtered),fft_filter_remap(sl),file_suffix(file_suffix) {}
    
    fft_filterer_t(fft_filterer_t &&)=default;
  };
  
  EXTERN_PROP std::vector<fft_filterer_t> fft_filterer;
  
  void start_hit(int ihit,bool skip=false);
}

#undef INIT_TO

#endif
