#ifndef _NEW_TYPES_DEFINITIONS_HPP
#define _NEW_TYPES_DEFINITIONS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#if HIGH_PREC == GMP_HIGH_PREC
 #include <gmpxx.h>
#endif

#ifdef USE_MPI
 #include <mpi.h>
#endif

#include <functional>
#include <math.h>
#include <stdint.h>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <sys/types.h>

#include "base/macros.hpp"

#ifdef SPI
 #include <spi/include/kernel/MU.h>
#endif

namespace nissa
{
  ///////////////// New types ///////////////////
  
  typedef int coords[NDIM];
  typedef double momentum_t[NDIM];
  
  typedef int intpair[2];
  typedef uint32_t checksum[2];
  
  typedef double complex[2];
  typedef complex quad_u1[NDIM];
  
  typedef complex spin[4];
  typedef complex color[3];
  typedef complex halfspin[2];
  
  typedef float single_complex[2];
  typedef single_complex single_color[NCOL];
  typedef single_color single_su3[NCOL];
  typedef single_color single_halfspincolor[2];
  typedef single_color single_spincolor[4];
  typedef single_su3 single_quad_su3[4];
  
  typedef spin colorspin[NCOL];
  typedef color spincolor[4];
  typedef color halfspincolor[2];
  
  typedef colorspin spincolorspin[4];
  typedef spincolorspin colorspincolorspin[NCOL];
  
  typedef spin spinspin[4];
  typedef spinspin colorspinspin[NCOL];
  
  typedef color su3[NCOL];
  typedef su3 quad_su3[NDIM];
  typedef su3 oct_su3[2*NDIM];
  
  typedef complex color2[2];
  typedef color2 su2[2];
  
  typedef colorspinspin su3spinspin[NCOL];
  
  typedef complex as2t[NDIM*(NDIM+1)/2];
  typedef su3 as2t_su3[NDIM*(NDIM+1)/2];
  typedef su3 opt_as2t_su3[4];
  
  typedef su3 squared_staples_t[NDIM][NDIM*(NDIM+1)/2];
  
  typedef squared_staples_t rectangular_staples_t;
  
#ifdef BGQ
  typedef complex bi_complex[2];
  typedef bi_complex bi_color[NCOL];
  typedef bi_color bi_su3[NCOL];
  typedef bi_su3 bi_oct_su3[8];
  typedef bi_color bi_spincolor[4];
  typedef bi_color bi_halfspincolor[2];
  typedef bi_complex bi_halfspin[2];
  typedef bi_su3 bi_opt_as2t_su3[4];
  
  typedef single_complex bi_single_complex[2];
  typedef bi_single_complex bi_single_color[NCOL];
  typedef bi_single_color bi_single_su3[NCOL];
  typedef bi_single_color bi_single_halfspincolor[2];
  typedef bi_single_color bi_single_spincolor[4];
  typedef bi_single_su3 bi_single_oct_su3[8];

#ifdef BGQ_EMU
  typedef bi_complex reg_bi_complex;
#else
  typedef vector4double reg_bi_complex;
#endif //BGQ_EMU
#endif //BGQ
  
#ifdef USE_MPI
#ifdef USE_MPI_IO
  typedef MPI_Offset ILDG_Offset;
  typedef MPI_File ILDG_File;
#else
  typedef off_t ILDG_Offset;
  typedef FILE* ILDG_File;
#endif
#endif
  
  //this is just for avoid misleading, but is nothing more that a spinspin
  typedef complex spin1field[4];
  typedef spin1field spin1prop[4];
  
  //quadruple precision float
  typedef double float_128[2];
  typedef float_128 complex_128[2];
  typedef complex_128 color_128[NCOL];
  typedef color_128 spincolor_128[4];

#ifdef BGQ
  typedef complex_128 bi_complex_128[2];
  typedef bi_complex_128 bi_color_128[NCOL];
  typedef bi_color_128 bi_su3_128[NCOL];
  typedef bi_su3_128 bi_oct_su3_128[8];
  typedef bi_color_128 bi_spincolor_128[4];
#endif

  //octpuple
  typedef double float_256[4];
  typedef double float_256_unr[5];
  struct float_256_class;
#if HIGH_PREC == GMP_HIGH_PREC
  typedef mpf_class float_high_prec_t;
#elif HIGH_PREC == NATIVE_HIGH_PREC
  typedef float_256_class float_high_prec_t;
#else
 #error Unknwon high_prec: HIGH_PREC
#endif

  //Random types
  enum rnd_t{RND_ALL_PLUS_ONE,RND_ALL_MINUS_ONE,RND_UNIF,RND_Z2,RND_Z4,RND_GAUSS};
  //Source type
  enum source_t{POINT_SOURCE,UNDILUTED_SOURCE,COLOR_DILUTED_SOURCE,SPIN_DILUTED_SOURCE,SPINCOLOR_DILUTED_SOURCE};
  //The three possibilities of quark computation
  enum hmc_force_piece{GAUGE_FORCE_ONLY,QUARK_FORCE_ONLY,BOTH_FORCE_PIECES};
  enum multistep_level{MACRO_STEP,MICRO_STEP};
  //Starting condition for a gauge conf
  enum start_conf_cond_t{UNSPEC_START_COND,HOT_START_COND,COLD_START_COND};
  //Boundary conditions
  enum boundary_cond_t{UNSPEC_BOUNDARY_COND,PERIODIC_BOUNDARY_COND,OPEN_BOUNDARY_COND};
  //Gauge action
  enum gauge_action_name_t{UNSPEC_GAUGE_ACTION,WILSON_GAUGE_ACTION,TLSYM_GAUGE_ACTION,IWASAKI_GAUGE_ACTION};
  //Basis for twisted mass
  enum tm_basis_t{WILSON_BASE,MAX_TWIST_BASE};
  
  typedef uint8_t coords_5D[5];
  
  ///////////////////// New structures ////////////////////
  
  //The structure for the random generator
  struct rnd_gen
  {
    int idum;
    int idum2;
    int iv[RAN2_NTAB];
    int iy;
  };
  
  //The structure for gamma matrix
  struct dirac_matr
  {
    int pos[4];
    complex entr[4];
  };
  
  //nissa vector
  struct nissa_vect
  {
    int64_t nel;
    int64_t size_per_el;
    
    char tag[NISSA_VECT_STRING_LENGTH];
    char type[NISSA_VECT_STRING_LENGTH];
    
    char file[NISSA_VECT_STRING_LENGTH];
    int line;
    
    nissa_vect *prev;
    nissa_vect *next;
    
    uint32_t flag;
    
    //padding to keep memory alignment
    char pad[(NISSA_VECT_ALIGNMENT-(2*sizeof(int64_t)+3*NISSA_VECT_STRING_LENGTH+sizeof(int)+2*sizeof(nissa_vect*)+sizeof(uint32_t))%NISSA_VECT_ALIGNMENT)%
	      NISSA_VECT_ALIGNMENT];
  };
}

#include "base/vectors.hpp"

namespace nissa
{
  struct su3_path
  {
    int ivol;
    int movements;
    su3 data;
    su3_path* next;
  };
  
  #if NCOL == 3
  //used to exponentiate for stouting
  struct hermitian_exp_ingredients
  {
    int sign;
    double c0,c1;
    double cu,su;
    double c2u,s2u;
    su3 Q,Q2;
    double u,w,theta;
    double xi0w;
    double cw;
    complex f[3];
  };
  #endif
  
  //ILDG header
  struct ILDG_header
  {
    uint32_t magic_no;
    uint16_t version;
    uint16_t mbme_flag;
    uint64_t data_length;
    char type[128];
  };
  
  //store messages
  struct ILDG_message
  {
    bool is_last;
    char *data;
    char *name;
    uint64_t data_length;
    ILDG_message *next;
  };
  
  //ILDG file view
  struct ILDG_File_view
  {
#ifdef USE_MPI
#ifdef USE_MPI_IO
    MPI_Datatype etype;
    MPI_Datatype ftype;
    MPI_Offset view_pos;
    MPI_Offset pos;
#endif
#endif
    char format[100];
  };
  
  ////////////////////// optimized 2 pts //////////////////////
  
  //match of forward and backward components id
  struct forw_back_comp_id_t: std::pair<uint8_t,uint8_t>
    {
      uint8_t &forw() {return this->first;}
      uint8_t &back() {return this->second;}
      forw_back_comp_id_t(uint8_t forw_comp_id,uint8_t back_comp_id)
	{
	  forw()=forw_comp_id;
	  back()=back_comp_id;
	}
    };
  
  //list of correlation function to which the combo contributes, and weight
  struct corr_id_weight_t: std::map<uint16_t,double> {};
  
  //Structure useful to compute a list of correlation function of the type
  // [Re/Im]Tr[G_sink*S_forw*G_sour*G5*S_back^+*G5]
  //we order S as sink-source-re/im.
  //To compute each correlation function we need to summ 32 contribution,
  //relative to the product of the real and imaginary components.
  //Two main patterns are possible, according to the the fact that we must 
  //compute the real or imaginary part of the correlation function,
  //and if one of the two dirac matrices are imaginary.
  //We compact the product of the various correlation functions, so that
  //we can compute only once the repeated products.
  struct two_pts_comp_t: std::map<forw_back_comp_id_t,corr_id_weight_t>
    {
      int ncorr;
      void add_sink_source_corr(uint16_t corr_id,double weight,int re_im,uint8_t sink_igamma,uint8_t sour_igamma);
      void print(FILE *fout=stdout);
      void summ_the_loc_forw_back_contractions(double *out,double *S_forw,double *S_back,int nel,int twall);
      void scream();
      std::map<int,std::string> corr_name;
      std::vector<int> pattern_list;
      
      void print_correlations_to_file(FILE *fout,double *corr);
    };
  
  //////////////////////////////////////////////////////////////////////
  
  // Type to hold the position of linear algebra output data (see "two_stage_computations" doc for explenations)
  struct two_stage_computation_pos_t
  {
    int *inter_fr_in_pos; //offset for intermediate result
    int *final_fr_inter_pos; //offset for final result from intermediate
    int *inter_fr_recv_pos; //offset for intermediate from nissa_recv_buf
    two_stage_computation_pos_t()
    {
      inter_fr_in_pos=final_fr_inter_pos=inter_fr_recv_pos=NULL;
    }
    void free()
    {
      if(inter_fr_in_pos!=NULL) nissa_free(inter_fr_in_pos);
      if(final_fr_inter_pos!=NULL) nissa_free(final_fr_inter_pos);
      if(inter_fr_recv_pos!=NULL) nissa_free(inter_fr_recv_pos);
    }
    ~two_stage_computation_pos_t() {free();}
  };
  
  //structure to store data for stouting
  struct stout_link_staples
  {
    su3 C;
    su3 Omega;
    su3 Q;
  };
  
  //parameters to stout
  struct stout_pars_t
  {
    int nlevels;
    double rho;
    
    int def_nlevels(){return 0;}
    double def_rho(){return 0;}
    
    int master_fprintf(FILE *fout,bool full=false);
    
    stout_pars_t() : nlevels(def_nlevels()),rho(def_rho()) {}
  };
  
  //parameters to ape smear
  struct ape_pars_t
  {
    int nlev;
    double alpha;
    
    ape_pars_t() : nlev(20),alpha(0.5) {}
  };
  
  //structure to cool
  struct cool_pars_t
  {
    gauge_action_name_t gauge_action;
    int nsteps;
    int overrelax_flag;
    double overrelax_exp;
    
    cool_pars_t() : gauge_action(WILSON_GAUGE_ACTION),nsteps(120),overrelax_flag(0),overrelax_exp(0.0) {}
  };
  
  //structure to adataptive stout
  struct adaptative_stout_pars_t
  {
    int nlevels;
    std::vector<double> rho;
    
    adaptative_stout_pars_t() : nlevels(0) {}
  };
  
  //structure to Wilson flow
  struct Wflow_pars_t
  {
    double T,dt;
    
    Wflow_pars_t() : T(10),dt(0.2) {}
  };
  
  //parameters to smooth a configuration
  struct smooth_pars_t
  {
    enum method_t{UNSPEC_SMOOTH_METHOD,COOLING,STOUTING,ADAPTATIVE_STOUTING,WFLOWING};
    
    //basic
    int flag;
    double meas_each;
    method_t method;
    //pars
    cool_pars_t cool_pars;
    stout_pars_t stout_pars;
    adaptative_stout_pars_t adaptative_stout_pars;
    Wflow_pars_t Wflow_pars;
    
    smooth_pars_t() : flag(false),method(COOLING) {}
  };
  
  //rational approximation
  struct rat_approx_t
  {
    char name[20];
    double minimum;
    double maximum;
    double maxerr;
    double cons;
    int num,den;
    int degree;
    double *poles;
    double *weights;
    void reset()
    {
      strcpy(name,"");
      minimum=maximum=cons=maxerr=0;
      degree=num=den=0;
      poles=weights=NULL;
    }
    rat_approx_t(){reset();}
  };
  
  //quark content
  struct quark_content_t
  {
    std::string name;
    int deg;
    double mass;
    double re_pot;
    double im_pot;
    double charge;
    
    std::string def_name(){return "quark";}
    int def_deg(){return 1;}
    double def_mass(){return 0.1;}
    double def_re_pot(){return 0;}
    double def_im_pot(){return 0;}
    double def_charge(){return 0;}
    
    int master_fprintf(FILE *fout,bool full=false);
    
    quark_content_t() : deg(def_deg()),mass(def_mass()),re_pot(def_re_pot()),im_pot(def_im_pot()),charge(def_charge()) {}
  };
  
  //parameters to smear a gauge conf when measuring gauge observables
  struct gauge_obs_temp_smear_pars_t
  {
    int use_hyp_or_ape_temp;
    
    //hyp pars
    double hyp_temp_alpha0;
    double hyp_temp_alpha1;
    double hyp_temp_alpha2;
    
    //ape temp
    double ape_temp_alpha;
    int nape_temp_iters;
    
    gauge_obs_temp_smear_pars_t() : use_hyp_or_ape_temp(0),hyp_temp_alpha0(1),hyp_temp_alpha1(1),hyp_temp_alpha2(0.5),ape_temp_alpha(0.5),nape_temp_iters(20) {}
  };
  
  //pars to compute polyakov loop
  struct poly_corr_meas_pars_t
  {
    int flag;
    std::string path;
    gauge_obs_temp_smear_pars_t gauge_smear_pars;
    int dir;
    
    int def_flag(){return 0;}
    std::string def_path(){return "luppoli";}
    int def_dir(){return 0;}
    
    int master_fprintf(FILE *fout,bool full=false);
    
    poly_corr_meas_pars_t() : flag(def_flag()),path(def_path()),dir(def_dir()) {}
  };
  
  //parameters to compute gauge observabls
  struct gauge_obs_meas_pars_t
  {
    int flag;
    std::string path;
    
    int def_flag(){return 0;}
    std::string def_path(){return "gauge_obs";}
    
    int master_fprintf(FILE *fout,bool full=false);
    
    gauge_obs_meas_pars_t() : flag(def_flag()),path(def_path()) {}
  };
  
  //parameters to compute the fermionic gran-mix
  struct fermionic_putpourri_meas_pars_t
  {
    int flag;
    std::string path;
    double residue;
    int compute_susc;
    int ncopies;
    int nhits;
    
    int def_flag(){return 0;}
    std::string def_path(){return "lavanda";}
    double def_residue(){return 1e-12;}
    int def_compute_susc(){return 0;}
    int def_ncopies(){return 1;}
    int def_nhits(){return 1;}
    
    int master_fprintf(FILE *fout,bool full=false);
    
    fermionic_putpourri_meas_pars_t() : flag(def_flag()),path(def_path()),residue(def_residue()),
					compute_susc(def_compute_susc()),ncopies(def_ncopies()),nhits(def_nhits()) {}
  };
  
  //parameters to compute the quark density and its high-order suscetivities
  struct quark_rendens_meas_pars_t
  {
    int flag;
    std::string path;
    double residue;
    int ncopies;
    int nhits;
    
    int def_flag(){return 0;}
    std::string def_path(){return "rende";}
    double def_residue(){return 1e-12;}
    int def_ncopies(){return 1;}
    int def_nhits(){return 1;}
    
    int master_fprintf(FILE *fout,bool full=false);
    
    quark_rendens_meas_pars_t() : flag(def_flag()),path(def_path()),residue(def_residue()),ncopies(def_ncopies()),nhits(def_nhits()) {}
  };
  
  //parameters to compute spin polarization
  struct spinpol_meas_pars_t
  {
    int flag;
    std::string path;
    double residue;
    int dir;
    int nhits;
    int use_ferm_conf_for_gluons;
    smooth_pars_t smooth_pars;
    
    int def_flag(){return false;}
    std::string def_path(){return "pollo";}
    double def_residue(){return 1e-12;}
    int def_dir(){return 1;}
    int def_nhits(){return 1;}
    int def_use_ferm_conf_for_gluons(){return 0;}
    
    spinpol_meas_pars_t() : flag(def_flag()),path(def_path()),residue(def_residue()),dir(def_dir()),nhits(def_nhits()),
			    use_ferm_conf_for_gluons(def_use_ferm_conf_for_gluons()) {}
  };
  
  //parameters to compute the magnetization
  struct magnetization_meas_pars_t
  {
    int flag;
    std::string path;
    double residue;
    int ncopies;
    int nhits;
    
    int def_flag(){return false;}
    std::string def_path(){return "magnetization";}
    double def_residue(){return 1e-12;}
    int def_ncopies(){return 1;}
    int def_nhits(){return 1;}
    
    int master_fprintf(FILE *fout,bool full=false);
    
    magnetization_meas_pars_t() : flag(def_flag()),path(def_path()),residue(def_residue()),ncopies(def_ncopies()),nhits(def_nhits()) {}
  };
  
  //parameters to compute time pseduoscalar correlator
  struct pseudo_corr_meas_pars_t
  {
    int flag;
    std::string path;
    double residue;
    int nhits;
    
    int def_flag(){return false;}
    std::string def_path(){return "pseudo_corr";}
    double def_residue(){return 1e-12;}
    int def_nhits(){return 1;}
    
    int master_fprintf(FILE *fout,bool full=false);
    
    pseudo_corr_meas_pars_t() : flag(def_flag()),path(def_path()),residue(def_residue()),nhits(def_nhits()) {}
  };
  
  //parameters to measure topology properties
  struct top_meas_pars_t
  {
    int flag;
    std::string path;
    smooth_pars_t smooth_pars;
    
    top_meas_pars_t() : flag(0),path("topocharge") {}
  };
  
  //holds temporal-spatial gauge smearing parameters
  struct gauge_obs_temp_spat_smear_pars_t
  {
    //hyp or ape in T direction and pars
    gauge_obs_temp_smear_pars_t gauge_temp_smear_pars;
    //ape spat
    double ape_spat_alpha;
    int nape_spat_levls,*nape_spat_iters;
    
    gauge_obs_temp_spat_smear_pars_t() : ape_spat_alpha(0.0),nape_spat_levls(0) {}
  };
  
  //parameters to measure all rectangles path
  struct all_rect_meas_pars_t
  {
    int flag;
    std::string path;
    gauge_obs_temp_spat_smear_pars_t smear_pars;
    int Tmin,Tmax,Dmin,Dmax;
    
    int def_flag(){return 0;}
    std::string def_path(){return "rectangles";}
    int def_Tmin(){return 3;}
    int def_Tmax(){return 9;}
    int def_Dmin(){return 1;}
    int def_Dmax(){return 9;}
     
    all_rect_meas_pars_t() : flag(def_flag()),path(def_path()),Tmin(def_Tmin()),Tmax(def_Tmax()),Dmin(def_Dmin()),Dmax(def_Dmax()) {}
  };
  
  //parameters to measure flux tube
  struct watusso_meas_pars_t
  {
    int flag;
    std::string path;
    gauge_obs_temp_spat_smear_pars_t smear_pars;
    int size_min,size_max,size_step,dmax;
    
    int def_flag(){return 0;}
    std::string def_path(){return "watusso";}
    int def_size_min(){return 7;}
    int def_size_max(){return 7;}
    int def_size_step(){return 1;}
    int def_dmax(){return 10;}
    
    watusso_meas_pars_t() : flag(def_flag()),path(def_path()),size_min(def_size_min()),size_max(def_size_max()),
			    size_step(def_size_step()),dmax(def_dmax()) {}
  };
  
  //storable vector
  ILDG_message* ILDG_string_message_append_to_last(ILDG_message *mess,const char *name,const char *data);
  template<class T> struct storable_vector_t : std::vector<T>
  {
    //append to last message
    ILDG_message *append_to_message_with_name(ILDG_message &mess,const char *name)
    {
      std::ostringstream os;
      os.precision(16);
      for(typename std::vector<T>::iterator it=this->begin();it!=this->end();it++) os<<*it<<" ";
      return ILDG_string_message_append_to_last(&mess,name,os.str().c_str());
    }
    //convert from a text message
    void convert_from_text(const char *data)
    {
      std::istringstream is(data);
      T temp;
      while(is>>temp) this->push_back(temp);
    }
    void convert_from_message(ILDG_message &mess)
    {convert_from_text(mess.data);}
    
    //read it from file
    void read_from_ILDG_file(ILDG_File fin,const char *tag);
  };
}

#include "metadynamics.hpp"

namespace nissa
{
  //parameters to add topological potential
  struct topotential_pars_t : meta_pars_t
  {
    int flag;
    double theta;
    stout_pars_t stout_pars;
    //methods inside opearations/su3_paths/topological_charge.cpp
    void store_if_needed(quad_su3 **conf,int iconf);
    
    int master_fprintf(FILE *fout,bool full=false);
    
    topotential_pars_t() : meta_pars_t(),flag(0),theta(0.0){}
  };
  
  //parameters to em field
  struct em_field_pars_t
  {
    int flag;
    
    //basic
    double E[3];
    double B[3];
    
    int master_fprintf(FILE *fout,bool full=false);
    
    em_field_pars_t() : flag(0) {for(int i=0;i<3;i++) E[i]=B[i]=0;}
  };
  
  //theory content
  struct theory_pars_t
  {
    double beta;
    int nflavs;
    quad_u1 ***backfield;
    quark_content_t *quark_content;
    gauge_action_name_t gauge_action_name;
    topotential_pars_t topotential_pars;
    stout_pars_t stout_pars;
    em_field_pars_t em_field_pars;
    
    theory_pars_t() : beta(6.0),nflavs(0) {}
  };
  
  //evolution parameters for hybrid monte carlo
  struct hmc_evol_pars_t
  {
    int skip_mtest_ntraj;
    double traj_length;
    double pf_action_residue;
    double md_residue;
    int nmd_steps;
    int ngauge_substeps;
    int *npseudo_fs;
    rat_approx_t *rat_appr;
    
    int master_fprintf(FILE *fout,bool full=false);
    
    hmc_evol_pars_t() : skip_mtest_ntraj(30),traj_length(1.0),pf_action_residue(1e-12),md_residue(1e-6),nmd_steps(13),ngauge_substeps(5) {}
  };
  
  //results of a unitarity check
  struct unitarity_check_result_t
  {
    int nbroken_links;
    double average_diff;
    double max_diff;
    
    unitarity_check_result_t () : nbroken_links(0),average_diff(0.0),max_diff(0.0) {}
  };
  
  //parameters for pure gauge theory
  struct pure_gauge_evol_pars_t
  {
    //wheter to use or not hmc
    int use_hmc;
    //basic hmc pars
    double traj_length;
    int nmd_steps;
    //acceleration parameters
    int use_Facc;
    double kappa;
    double residue;
    //number of hb sweeps and hits per link
    int nhb_sweeps;
    int nhb_hits;
    //the same for overrelax
    int nov_sweeps;
    int nov_hits;
    
    pure_gauge_evol_pars_t() : use_hmc(0),traj_length(1.0),nmd_steps(13),use_Facc(0),kappa(0.0),residue(1e-12),nhb_sweeps(1),nhb_hits(1),nov_sweeps(3),nov_hits(3) {}
  };
  
  struct evol_pars_t
  {
    hmc_evol_pars_t hmc_evol_pars;
    pure_gauge_evol_pars_t pure_gauge_evol_pars;
  };
  
#ifdef USE_MPI
  //out and in buffer
  struct comm_t
  {
    //bgq specific structures, in alternative to ordinary MPI
#ifdef SPI
    //descriptors
    MUHWI_Descriptor_t *descriptors;
    MUHWI_Destination spi_dest[8];
#else
    //destinations and source ranks
    int send_rank[8],recv_rank[8];
#    //requests and message
    MPI_Request requests[16];
    int nrequest,imessage;
#endif
    
    //communication in progress
    int comm_in_prog;
    //local size
    uint64_t nbytes_per_site;
    //size of the message
    uint64_t tot_mess_size;
    //offsets
    int send_offset[8],message_length[8],recv_offset[8];
    
    //constructor
    bool initialized;
    comm_t(){initialized=false;}
  };
#endif
}
#endif
