#ifndef _PROP_HPP
#define _PROP_HPP

#include <filesystem>
#include <set>

#include "conf.hpp"
#include "meslep.hpp"
#include "pars.hpp"

#ifndef EXTERN_PROP
# define EXTERN_PROP extern
#define INIT_TO(A)
#else
# define INIT_TO(A) =A
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
    
    double kappa{};
    Momentum kappa_asymm{};
    double mass{};
    int r{};
    double charge{};
    Momentum theta{};
    
    insertion_t insertion{};
    std::vector<source_term_t> source_terms{};
    int tins{};
    double residue{};
    bool store{};
    
    char ext_field_path[33]{};
    
    rnd_t noise_type{};
    double ori_source_norm2{};
    
    //spincolor data
    LxField<spincolor>** sp{};
    
    CUDA_HOST_AND_DEVICE
    const LxField<spincolor>& operator[](const int& i) const
    {
      return *(sp[i]);
    }
    
    CUDA_HOST_AND_DEVICE
    LxField<spincolor>& operator[](const int& i)
    {
      return *(sp[i]);
    }
    
    void alloc_storage()
    {
      if(sp==nullptr)
	{
	  sp=new LxField<spincolor>*[nso_spi*nso_col];
	  for(int i=0;i<nso_spi*nso_col;i++)
	    sp[i]=new LxField<spincolor>("sp",WITH_HALO);
	}
    }
    
    //initialize as a propagator
    void init_as_propagator(const insertion_t& _insertion,
			    const std::vector<source_term_t>& _source_terms,
			    const int& _tins,
			    const double& _residue,
			    const double& _kappa,
			    const double* _kappa_asymm,
			    const double& _mass,
			    const char *_ext_field_path,
			    const int& _r,
			    const double& _charge,
			    const Momentum& _theta,
			    const bool& _store)
    {
      is_source=false;
      
      kappa=_kappa;
      for(int mu=0;mu<NDIM;mu++)
	kappa_asymm[mu]=_kappa_asymm[mu];
      
      mass=_mass;
      
      r=_r;
      
      charge=_charge;
      
      for(int mu=0;mu<NDIM;mu++)
	theta[mu]=_theta[mu];
      
      insertion=_insertion;
      
      source_terms=_source_terms;
      
      tins=_tins;
      
      strncpy(ext_field_path,_ext_field_path,32);
      
      residue=_residue;
      
      store=_store;
      
      if(is_photon_ins(insertion))
	need_photon=true;
      
      alloc_storage();
    }
    
    //initialize as a source
    void init_as_source(const rnd_t& _noise_type,
			const int& _tins,
			const int& _r,
			const bool& _store)
    {
      is_source=true;
      
      noise_type=_noise_type;
      
      tins=_tins;
      
      r=_r;
      
      store=_store;
      
      alloc_storage();
    }
    
    void free_storage()
    {
      if(sp)
	{
	  for(int i=0;i<nso_spi*nso_col;i++)
	    delete sp[i];
	  delete[] sp;
	  
	  sp=nullptr;
	}
    }
    
    qprop_t(const insertion_t& insertion,
	    const std::vector<source_term_t>& source_terms,
	    const int& tins,
	    const double& residue,
	    const double& kappa,
	    const double* kappa_asymm,
	    const double& mass,
	    const char* ext_field_path,
	    const int& r,
	    const double& charge,
	    const Momentum& theta,
	    const bool& store)
    {
      init_as_propagator(insertion,source_terms,tins,residue,kappa,kappa_asymm,mass,ext_field_path,r,charge,theta,store);
    }
    
    qprop_t(const rnd_t& noise_type,
	    const int& tins,
	    const int& r,
	    const bool& store)
    {
      init_as_source(noise_type,tins,r,store);
    }
    
    qprop_t()
    {
      is_source=0;
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
  
  EXTERN_PROP LxField<spincolor> *loop_source;
  
  inline void allocate_loop_source()
  {
    loop_source=new LxField<spincolor>("loop_source",WITH_HALO);
  }
  
  inline void free_loop_source()
  {
    delete loop_source;
  }
  
  EXTERN_PROP int nsource_tot INIT_TO(0),nphoton_prop_tot INIT_TO(0);
  EXTERN_PROP double source_time INIT_TO(0),photon_prop_time INIT_TO(0),lepton_prop_time INIT_TO(0);
  
  EXTERN_PROP int load_photons INIT_TO(false);
  EXTERN_PROP int store_photons INIT_TO(false);
  
  CUDA_MANAGED EXTERN_PROP LxField<spin1field> *photon_field;
  
  EXTERN_PROP LxField<spin1field> *photon_phi;
  
  EXTERN_PROP LxField<spin1field> *photon_eta;
  
  void allocate_photon_fields();
  void free_photon_fields();
  CUDA_MANAGED EXTERN_PROP spinspin *temp_lep;
  
  void generate_original_sources(const int& ihit,
				 const bool& skipOnly);
  
  void generate_original_source(qprop_t& sou,
				const bool& skipOnly);
  
  void generate_quark_propagator(std::string& name,qprop_t& q,int ihit);
  void generate_photon_source(LxField<spin1field>& photon_eta);
  
  void generate_source(insertion_t inser,int r,double charge,double kappa,const Momentum& theta,spincolor *ori,int t);
  
  void generate_quark_propagators(const int& ihit);
  
  void generate_photon_stochastic_propagator(const int& ihit);
  
  //CUDA_HOST_AND_DEVICE void get_antineutrino_source_phase_factor(complex out,const int ivol,const int ilepton,const Momentum bc);
  void generate_lepton_propagators();
  void propagators_fft(const int& ihit);
  
  /// Multiply the configuration for an additional u(1) field, defined as exp(-i e q A /3)
  void add_photon_field_to_conf(LxField<quad_su3>& conf,
				const double& charge);
  
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
  
  inline void start_hit(int ihit,bool skip=false)
  {
    MASTER_PRINTF("\n=== Hit %d/%d ====\n",ihit+1,nhits);
    if(use_new_generator)
      {
	for(int mu=0;mu<NDIM;mu++)
	  {
	    using C=double[1];
	    C c;
	    field_rng_stream.drawScalar(c);
	    source_coord[mu]=c[0]*glbSize[mu];
	  }
      }
    else
      source_coord=generate_random_coord();
    
    if(stoch_source) MASTER_PRINTF(" source time: %d\n",source_coord[0]);
    else             MASTER_PRINTF(" point source coords: %d %d %d %d\n",source_coord[0],source_coord[1],source_coord[2],source_coord[3]);
    if(need_photon)
      {
	if(skip)
	  generate_photon_source(*photon_eta);
	else
	  generate_photon_stochastic_propagator(ihit);
      }
    generate_original_sources(ihit,skip);
  }
  
  inline void generate_propagators(int ihit)
  {
    if(nquark_lep_combos) generate_lepton_propagators();
    generate_quark_propagators(ihit);
  }
  
  template <typename T>
  struct ReadWriteRealVector
  {
    LxField<T>& v;
    
    std::string path;
    
    FILE* fastFile;
    
    ReadWriteRealVector(LxField<T>& v,
			const std::string& _path) :
      v(v),
      path(_path)
    {
      if(fast_read_write_vectors)
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
	CRASH("Unable to open path %s with mode %s",path.c_str(),mode);
    }
    
    void cleanFiles()
    {
      std::filesystem::remove(path);
    }
    
    void fastRead()
    {
      CRASH("reimplement");
      fastOpen("r");
      
      // if(fread(v,sizeof(T),locVol,fastFile)!=(size_t)locVol)
      // 	CRASH("Problem reading %s",path.c_str());
      
      fclose(fastFile);
    }
    
    void read()
    {
      MASTER_PRINTF("Reading %s\n",path.c_str());
      
      if(fast_read_write_vectors)
	{
	  CRASH("reimplement");
	  fastOpen("r");
	  
	  // if(fread(v,sizeof(T),locVol,fastFile)!=locVol)
	  //   CRASH("Problem reading %s",path.c_str());
	  
	  fclose(fastFile);
	}
      else
	read_real_vector(v,path,"scidac-binary-data");
    }
    
    void fastWrite()
    {
	  CRASH("reimplement");
      fastOpen("w");
      
      // if(fwrite(v,sizeof(T),locVol,fastFile)!=(size_t)locVol)
      // 	CRASH("Problem writing %s",path.c_str());
      
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
// 	CRASH("Problem writing %s, total written: %zu when %zu expected, seek worked: %d",path.c_str(),written,totData,(int)seekingError);
      
      fclose(fastFile);
    }
    
    void write()
    {
      MASTER_PRINTF("Writing %s, %zu %p\n",path.c_str(),sizeof(T),&v);
      
      if(fast_read_write_vectors)
	{
	  fastOpen("w");
	  if(std::remove_reference_t<decltype(v)>::fieldLayout!=FieldLayout::CPU)
	    CRASH("not supported");
	  //hack
	  if(fwrite(v._data,sizeof(T),locVol,fastFile)!=(size_t)locVol)
	    CRASH("Problem writing %s",path.c_str());
	  
	  fclose(fastFile);
	}
      else
	write_real_vector(path,v,"scidac-binary-data");
    }
  };
}

#undef INIT_TO

#endif
