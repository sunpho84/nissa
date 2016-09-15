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

#define QUARK_BOUND_COND 1

namespace nissa
{
  EXTERN_PARS int diluted_spi_source,diluted_col_source;
  EXTERN_PARS int nso_spi,nso_col;
  
  //convention on gospel
  const int follow_chris=0,follow_nazario=1;
  
  //define types of quark propagator used
  const int nins_kind=8;
  enum insertion_t{                       ORI_SOU,  SCALAR,  PSEUDO,  PHOTON,  PHOTON_ETA,  PHOTON_PHI,  TADPOLE };//,  VECTOR};
  const insertion_t ins_list[nins_kind]={ ORI_SOU , SCALAR , PSEUDO , PHOTON , PHOTON_ETA , PHOTON_PHI , TADPOLE };//,  VECTOR };
  const char ins_name[nins_kind][20]=   {"ORI_SOU","SCALAR","PSEUDO","PHOTON","PHOTON_ETA","PHOTON_PHI","TADPOLE"};//, "VECTOR"};
  const char ins_tag[nins_kind]=        {'O'      ,'S'     ,'P'     ,'F'     ,'A'         ,'C'         ,'T'      };//, 'V'};
  inline insertion_t ins_from_tag(const char tag)
  {
    int i=0;
    while(i<nins_kind && ins_tag[i]!=tag) i++;
    if(i>=nins_kind) crash("unable to find tag %c",tag);
    return ins_list[i];
  }
  //sign of the lepton momentum
  const int norie=2;
  const int sign_orie[2]={-1,+1};
  
  EXTERN_PARS int nr_lep;
  const int nins=3;
  const int nlins=2;
  const int nrev=2;
  
  EXTERN_PARS int twisted_run,clover_run;
  EXTERN_PARS tm_basis_t base;
  EXTERN_PARS double glb_cSW;
  
  EXTERN_PARS tm_quark_info *leps;
  
  EXTERN_PARS gauge_info photon;
  EXTERN_PARS double tadpole[NDIM];
  
  EXTERN_PARS int nquark_lep_combos;
  EXTERN_PARS double *lep_energy,*neu_energy;
  
  void read_input_preamble();
  void read_mes2pts_contr_pars();
  void read_mes2pts_contr_gamma_list();
  void read_meslep_contr_pars();
  void read_bar2pts_contr_pars();
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
  
  //initialize the dilutions
  inline void read_set_dilutions_if_stoch_source(int stoch_source)
  {
    int dil_spin=1,dil_col=1;
    if(stoch_source)
      {
	read_str_int("DilutedSpin",&dil_spin);
	read_str_int("DilutedColor",&dil_col);
      }
    set_diluted_spin(dil_spin);
    set_diluted_color(dil_col);
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
  EXTERN_PARS int follow_chris_or_nazario INIT_TO(follow_nazario);
  inline void read_gospel_convention()
  {read_str_int("FollowChrisOrNazario",&follow_chris_or_nazario);}
  
  //perform a random gauge transformation
  EXTERN_PARS int rnd_gauge_transform INIT_TO(0);
  inline void read_random_gauge_transform()
  {read_str_int("RandomGaugeTransform",&rnd_gauge_transform);}
  
  //local pion or muon current?
  EXTERN_PARS int loc_hadr_curr INIT_TO(false);
  inline void read_loc_hadr_curr()
  {read_str_int("LocHadrCurr",&loc_hadr_curr);}
  EXTERN_PARS int loc_muon_curr INIT_TO(false);
  inline void read_loc_muon_curr()
  {read_str_int("LocMuonCurr",&loc_muon_curr);}
  
  //stochastic sources
  EXTERN_PARS int stoch_source INIT_TO(0);
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
  
  //gaussian smearing pars
  EXTERN_PARS int gaussian_smearing_niters;
  EXTERN_PARS double gaussian_smearing_kappa;
  inline void read_gaussian_smearing_pars()
  {  
    read_str_double("GaussianSmearingKappa",&gaussian_smearing_kappa);
    read_str_int("GaussianSmearingNiters",&gaussian_smearing_niters);
  }
  
  //read all smearing pars
  inline void read_smearing_pars()
  {
    read_ape_smearing_pars();
    read_gaussian_smearing_pars();
  }
}

#undef EXTERN_PARS
#undef INIT_TO

#endif
