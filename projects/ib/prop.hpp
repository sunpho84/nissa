#ifndef _PROP_HPP
#define _PROP_HPP

#include <filesystem>
#include <set>

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
    bool is_source{false};
    
    double kappa{0};
    double kappa_asymm[NDIM]{};
    double mass{0};
    int r{0};
    double charge{0};
    momentum_t theta{0};
    
    insertion_t insertion{};
    std::vector<source_term_t> source_terms{};
    int tins{};
    double residue{};
    bool store{};
    
    char ext_field_path[33]{};
    
    rnd_t noise_type{};
    double ori_source_norm2{};
    
    //spincolor data
    spincolor** sp{};
    
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
    
    void alloc_storage()
    {
      if(sp==nullptr)
	{
	  sp=nissa_malloc("sp",nso_spi*nso_col,spincolor*);
	  for(int i=0;i<nso_spi*nso_col;i++)
	    sp[i]=nissa_malloc("sp",locVol+bord_vol,spincolor);
	}
    }
    
    //initialize as a propagator
    void init_as_propagator(insertion_t _insertion,const std::vector<source_term_t>& _source_terms,int _tins,double _residue,double _kappa,const double* _kappa_asymm,double _mass,char *_ext_field_path,int _r,double _charge,const momentum_t& _theta,bool _store)
    {
      is_source=false;
      
      kappa=_kappa;
      for(int mu=0;mu<NDIM;mu++) kappa_asymm[mu]=_kappa_asymm[mu];
      mass=_mass;
      r=_r;
      charge=_charge;
      for(int mu=0;mu<NDIM;mu++) theta[mu]=_theta[mu];
      insertion=_insertion;
      source_terms=_source_terms;
      tins=_tins;
      strncpy(ext_field_path,_ext_field_path,32);
      residue=_residue;
      store=_store;
      
      if(is_photon_ins(insertion)) need_photon=true;
    }
    
    //initialize as a source
    void init_as_source(rnd_t _noise_type,int _tins,int _r,bool _store)
    {
      is_source=true;
      
      noise_type=_noise_type;
      tins=_tins;
      r=_r;
      store=_store;
    }
    
    void free_storage()
    {
      if(sp)
	{
	  for(int i=0;i<nso_spi*nso_col;i++)
	    nissa_free(sp[i]);
	  nissa_free(sp);
	}
    }
    
    ~qprop_t()
    {
      free_storage();
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
  EXTERN_PROP std::set<std::string> propsNeededToContr;
  EXTERN_PROP spinspin **L;
  
  EXTERN_PROP int nlprop;
  
  EXTERN_PROP std::vector<std::string> ori_source_name_list;
  EXTERN_PROP spincolor *loop_source;
  inline void allocate_loop_source(){loop_source=nissa_malloc("loop_source",locVol+bord_vol,spincolor);}
  inline void free_loop_source(){nissa_free(loop_source);}
  
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
  
  void get_qprop(spincolor *out,spincolor *in,double kappa,double mass,int r,double q,double residue,const momentum_t& theta);
  void insert_external_loc_source(spincolor *out,spin1field *curr,spincolor *in,int t,bool *dirs);
  void insert_external_source(spincolor *out,quad_su3 *conf,spin1field *curr,spincolor *ori,int t,int r,bool *dirs,int loc);
  void generate_quark_propagator(std::string& name,qprop_t& q,int ihit);
  void generate_photon_source(spin1field *photon_eta);
  void generate_source(insertion_t inser,int r,double charge,double kappa,const momentum_t& theta,spincolor *ori,int t);
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
  
  template <typename T>
  struct ReadWriteRealVector
  {
    T* v;
    std::string path;
    bool fast;
    FILE* fastFile;
    
    ReadWriteRealVector(T* v,const std::string& _path,const bool fast=fast_read_write_vectors) :
      v(v),path(_path),fast(fast)
    {
      if(fast)
	path+="_rank"+std::to_string(rank);
    }
    
    bool canLoad() const
    {
      return file_exists(path);
    }
    
    void fastOpen(const char* mode)
    {
      fastFile=fopen(path.c_str(),mode);
      if(fastFile==nullptr)
	crash("Unable to open path %s with mode %s",path.c_str(),mode);
    }
    
    void cleanFiles()
    {
      std::filesystem::remove(path);
    }
    
    void fastRead()
    {
      fastOpen("r");
      
      if(fread(v,sizeof(T),locVol,fastFile)!=(size_t)locVol)
	crash("Problem reading %s",path.c_str());
      
      fclose(fastFile);
    }
    
    void read()
    {
      master_printf("Reading %s\n",path.c_str());
      
      if(fast_read_write_vectors)
	fastRead();
      else
	read_real_vector(v,path,"scidac-binary-data");
    }
    
    void fastWrite()
    {
      fastOpen("w");
      
      if(fwrite(v,sizeof(T),locVol,fastFile)!=(size_t)locVol)
	crash("Problem writing %s",path.c_str());
      
//       size_t written=0;
//       bool seekingError=false;
//       const size_t totData=locVol*sizeof(T);
// #pragma omp parallel reduction(+:written) reduction(||:seekingError)
//       {
// 	const size_t nThr=omp_get_num_threads();
// 	const size_t iThr=omp_get_thread_num();
// 	const size_t maxThrData=(totData+nThr-1)/nThr;
// 	const size_t begData=maxThrData*iThr;
// 	const size_t endData=std::min(begData+maxThrData,totData);
// 	const size_t thrData=endData-begData;
// 	seekingError|=fseek(fastFile,begData,SEEK_SET);
// 	written+=fwrite(v,1,thrData,fastFile);
//       }
//       if(written!=totData or seekingError)
// 	crash("Problem writing %s, total written: %zu when %zu expected, seek worked: %d",path.c_str(),written,totData,(int)seekingError);
      
      fclose(fastFile);
    }
    
    void write()
    {
      master_printf("Writing %s, %zu %p\n",path.c_str(),sizeof(T),v);
      
      if(fast_read_write_vectors)
	fastWrite();
      else
	write_real_vector(path,v,64,"scidac-binary-data");
    }
  };
}

#undef INIT_TO

#endif
