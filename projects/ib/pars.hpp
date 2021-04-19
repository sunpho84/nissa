#ifndef _PARS_HPP
#define _PARS_HPP

#include <nissa.hpp>

#include "conf.hpp"

#ifndef EXTERN_PARS
 #define EXTERN_PARS extern
 #define INIT_TO(VAL)
#else
 #define INIT_TO(VAL) =VAL
#endif

namespace nissa
{
#define PERIODIC_BC 0
#define ANTIPERIODIC_BC 1
  
  //Twisted run
  EXTERN_PARS int twisted_run;
  EXTERN_PARS tm_basis_t base;
  inline void read_twisted_run()
  {
    read_str_int("TwistedRun",&twisted_run);
    if(!twisted_run) base=WILSON_BASE;
    else             base=MAX_TWIST_BASE;
  }
  
  //Clover run
  EXTERN_PARS int clover_run;
  EXTERN_PARS double glb_cSW;
  inline void read_clover_run()
  {
    read_str_int("CloverRun",&clover_run);
    if(clover_run) read_str_double("cSW",&glb_cSW);
  }
  
  CUDA_MANAGED EXTERN_PARS double temporal_bc INIT_TO(ANTIPERIODIC_BC);
  
  CUDA_MANAGED EXTERN_PARS int diluted_spi_source,diluted_col_source,diluted_spat_source;
  CUDA_MANAGED EXTERN_PARS int nso_spi,nso_col;
  CUDA_MANAGED EXTERN_PARS GlbCoords source_coord;

  template <typename C>
  CUDA_HOST_DEVICE inline
  C rel_coord_of_glb_coord(const C& c,const Direction& mu)
  {
    return (glbSize(mu)+c-source_coord(mu))%glbSize(mu);
  }
  
  inline GlbCoord rel_time_of_glb_time(const GlbCoord& t)
  {
    return rel_coord_of_glb_coord(t,0);
  }
  
  CUDA_HOST_DEVICE inline GlbCoord rel_coord_of_loclx(const LocLxSite& loclx,const Direction& mu)
  {
    return rel_coord_of_glb_coord(glbCoordOfLoclx(loclx,mu),mu);
  }
  
  CUDA_HOST_DEVICE inline GlbCoord rel_time_of_loclx(const LocLxSite& loclx)
  {
    return rel_coord_of_loclx(loclx,timeDirection);
  }
  
  //convention on gospel
  const int follow_chris=0,follow_nazario=1;
  
  //define types of quark propagator used
  const int nins_kind=23;
  enum insertion_t{                       PROP , SCALAR , PSEUDO , PHOTON , PHOTON_ETA , PHOTON_PHI , TADPOLE , CVEC , CVEC0 , CVEC1 , CVEC2 , CVEC3 , PHOTON0 , PHOTON1 , PHOTON2 , PHOTON3 , WFLOW , BACK_WFLOW , SMEARING , ANYSM , PHASING , EXT_FIELD , GAMMA };
  const insertion_t ins_list[nins_kind]={ PROP , SCALAR , PSEUDO , PHOTON , PHOTON_ETA , PHOTON_PHI , TADPOLE , CVEC , CVEC0 , CVEC1 , CVEC2 , CVEC3 , PHOTON0 , PHOTON1 , PHOTON2 , PHOTON3 , WFLOW , BACK_WFLOW , SMEARING, ANYSM, PHASING , EXT_FIELD , GAMMA };
  const char ins_name[nins_kind][20]=   {"PROP","SCALAR","PSEUDO","PHOTON","PHOTON_ETA","PHOTON_PHI","TADPOLE","CVEC","CVEC0","CVEC1","CVEC2","CVEC3","PHOTON0","PHOTON1","PHOTON2","PHOTON3","WFLOW","BACK_WFLOW","SMEARING","ANYSM","PHASING","EXT_FIELD","GAMMA"};
  const char ins_tag[nins_kind][3]=    {"-"   ,"S"     ,"P"     ,"F"     ,"A"         ,"C"         ,"T"      ,"V"   ,"V0"   ,"V1"   ,"V2"   ,"V3"   ,"F0"     ,"F1"     ,"F2"     ,"F3"     ,"WF"   ,"BF"       ,"SM"     ,"AN"     ,"PH"     ,"X"    ,"G"    };
  inline insertion_t ins_from_tag(const char *tag)
  {
    int i=0;
    while(i<nins_kind and strcasecmp(ins_tag[i],tag)) i++;
    if(i>=nins_kind)
      {
	master_fprintf(stderr,"unable to find tag %s, use one in the list:\n",tag);
	for(i=0;i<nins_kind;i++) master_fprintf(stderr," %s\n",ins_tag[i]);
	crash("see previous error");
      }
    return ins_list[i];
  }
  
  inline int is_smearing_ins(insertion_t ins)
  {
    return
      (ins==SMEARING) or
      (ins==ANYSM);
  }
  
  inline int is_photon_ins(insertion_t ins)
  {
    return
      (ins==PHOTON) or
      (ins==PHOTON_ETA) or
      (ins==PHOTON_PHI);
  }
  
  EXTERN_PARS gauge_info photon;
  EXTERN_PARS Momentum tadpole;
  
  //holds the range of FFT moms
  struct fft_mom_range_t
  {
    GlbCoords offs;
    GlbCoords width;
    
    fft_mom_range_t()
    {
    }
    
    fft_mom_range_t(const fft_mom_range_t& oth)
    {
      offs.nastyCopy(oth.offs);
      width.nastyCopy(oth.width);
    }
  };
  
  //list of propagators to fft
  EXTERN_PARS std::vector<std::string> fft_prop_list;
  
  void read_input_preamble();
  void read_mes2pts_contr_pars();
  void read_mes2pts_contr_gamma_list();
  void read_bar2pts_contr_pars();
  void read_handcuffs_contr_pars();
  void read_fft_prop_pars();
  void read_photon_pars();
  
  //set or not diluted the spin
  inline void set_diluted_spin(int s)
  {
    diluted_spi_source=s;
    if(diluted_spi_source) nso_spi=NDIRAC;
    else nso_spi=1;
  }
  
  //set or not diluted the color
  inline void set_diluted_color(int c)
  {
    diluted_col_source=c;
    if(diluted_col_source) nso_col=NCOL;
    else nso_col=1;
  }
  
  //set the dilution in space
  inline void set_diluted_space(int s)
  {diluted_spat_source=s;}
  
  //initialize the dilutions
  inline void read_set_dilutions_if_stoch_source(int stoch_source)
  {
    int dil_spin=1,dil_col=1,dil_spa=1;
    if(stoch_source)
      {
	read_str_int("DilutedSpin",&dil_spin);
	read_str_int("DilutedColor",&dil_col);
	read_str_int("DilutedSpace",&dil_spa);
      }
    set_diluted_spin(dil_spin);
    set_diluted_color(dil_col);
    set_diluted_space(dil_spa);
  }
  
  //initialize the random generator with the read seed
  inline void read_seed_start_random()
  {
    int seed;
    read_str_int("Seed",&seed);
    start_loc_rnd_gen(seed);
  }
  
  //flag to simulate in the free theory
  EXTERN_PARS int free_theory INIT_TO(false);
  inline void read_free_theory_flag()
  {read_str_int("FreeTheory",&free_theory);}
  
  //flag to make the muon with or without the external line
  CUDA_MANAGED EXTERN_PARS int follow_chris_or_nazario INIT_TO(follow_nazario);
  inline void read_gospel_convention()
  {read_str_int("FollowChrisOrNazario",&follow_chris_or_nazario);}
  
  //perform a random gauge transformation
  EXTERN_PARS int rnd_gauge_transform INIT_TO(0);
  inline void read_random_gauge_transform()
  {read_str_int("RandomGaugeTransform",&rnd_gauge_transform);}
  
  //perform a Landau gauge fixing
  EXTERN_PARS int Landau_gauge_fix_flag INIT_TO(0);
  EXTERN_PARS LC_gauge_fixing_pars_t gauge_fixing_pars;
  inline void read_Landau_gauge_fix()
  {
    read_str_int("LandauGaugeFix",&Landau_gauge_fix_flag);
    if(Landau_gauge_fix_flag) read_LC_gauge_fixing_pars(gauge_fixing_pars);
  }
  
  //store the conf?
  EXTERN_PARS int store_conf INIT_TO(0);
  inline void read_store_conf()
  {read_str_int("StoreConf",&store_conf);}
  
  //local pion or muon current?
  EXTERN_PARS int loc_hadr_curr INIT_TO(false);
  inline void read_loc_hadr_curr()
  {
    read_str_int("LocHadrCurr",&loc_hadr_curr);
  }

  EXTERN_PARS int loc_muon_curr INIT_TO(false);
  inline void read_loc_muon_curr()
  {read_str_int("LocMuonCurr",&loc_muon_curr);}
  
  //stochastic sources
  CUDA_MANAGED EXTERN_PARS int stoch_source INIT_TO(0);
  inline void read_stoch_source()
  {read_str_int("StochSource",&stoch_source);}
  
  //number of hits
  EXTERN_PARS int nhits INIT_TO(1);
  inline void read_nhits()
  {read_str_int("NHits",&nhits);}
  
  //number of configurations
  inline void read_ngauge_conf()
  {read_str_int("NGaugeConf",&ngauge_conf);}
  
  //ape smearing pars
  EXTERN_PARS int ape_smearing_niters;
  EXTERN_PARS double ape_smearing_alpha;
  inline void read_ape_smearing_pars()
  {
    read_str_double("ApeSmearingAlpha",&ape_smearing_alpha);
    read_str_int("ApeSmearingNiters",&ape_smearing_niters);
  }
  
  //read if isothrope theta or not
  EXTERN_PARS bool iso_theta INIT_TO(false);
  inline void read_iso_theta()
  {
    char theta_tag[1024];
    read_str(theta_tag,1024);
    
    if(strcasecmp(theta_tag,"Theta")==0) iso_theta=true;
    else
      if(strcasecmp(theta_tag,"ThetaX")==0)
	{
	  iso_theta=false;
	  expect_str("ThetaY");
	  expect_str("ThetaZ");
	}
      else crash("Unknown theta tag: %s",theta_tag);
  }
  
  //handle to stop
  EXTERN_PARS std::string stop_path INIT_TO("stop");
  
  //read the theta, iso or not
  inline void read_theta(Momentum& theta)
  {
    if(iso_theta)
      {
	read_double(&theta(xDirection));
	master_printf("Read variable 'Theta' with value: %lg\n",theta(xDirection));
	for(Direction mu=2;mu<NDIM;mu++)
	  theta(mu)=theta(xDirection);
      }
    else
      FOR_ALL_SPATIAL_DIRECTIONS(mu)
	{
	  read_double(&theta(mu));
	  master_printf("Read variable 'Theta[%d]' with value: %lg\n",mu,theta(mu));
	}
  }
  
  //lock the file to ensure single disk access
  EXTERN_PARS lock_file_t<uint64_t> lock_file;
}

#undef EXTERN_PARS
#undef INIT_TO

#endif
