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
  //for baryons
  const int nsm_sink=2;
  
  //convention on gospel
  const int follow_chris=0,follow_nazario=1;
  
  //define types of quark propagator used
  const int nins_kind=7;
  enum insertion_t{                    ORIGINAL,  SCALAR,  PSEUDO,  PHOTON,  PHOTON_ETA,  PHOTON_PHI,  TADPOLE};//,  VECTOR};
  const char ins_name[nins_kind][20]={"ORIGINAL","SCALAR","PSEUDO","PHOTON","PHOTON_ETA","PHOTON_PHI","TADPOLE"};//, "VECTOR"};
  
  //sign of the lepton momentum
  const int norie=2;
  const int sign_orie[2]={-1,+1};
  
  const int nins=3;
  const int nlins=2;
  const int nrev=2;
  
  EXTERN_PARS int pure_wilson;
  EXTERN_PARS tm_basis_t base;
  EXTERN_PARS double kappa;
  
  EXTERN_PARS tm_quark_info *leps;
  
  EXTERN_PARS int nqmass,nr;
  EXTERN_PARS double *qmass,*qkappa,*residue;
  
  EXTERN_PARS gauge_info photon;
  EXTERN_PARS double tadpole[NDIM];
  
  EXTERN_PARS int nleptons;
  EXTERN_PARS double *lep_energy,*neu_energy;
  
  void read_input_preamble();
  void read_photon_pars();
  
  //initialize the random generator with the read seed
  inline void read_seed_start_random()
  {
    int seed;
    read_str_int("Seed",&seed);
    start_loc_rnd_gen(seed);
  }
  
  //flag to store propagator0
  EXTERN_PARS int store_prop0 INIT_TO(false);
  inline void read_store_prop0_flag()
  {read_str_int("StoreProp0",&store_prop0);}
  
  //flag to simulate in the free theory
  EXTERN_PARS int free_theory INIT_TO(false);
  inline void read_free_theory_flag()
  {read_str_int("FreeTheory",&free_theory);}
  
  //flag to make the muon with or without the external line
  EXTERN_PARS int follow_chris_or_nazario INIT_TO(follow_nazario);
  inline void read_gospel_convention()
  {read_str_int("FollowChrisOrNazario",&follow_chris_or_nazario);}
  
  //noise type
  EXTERN_PARS int noise_type;
  inline void read_noise_type()
  {read_str_int("NoiseType",&noise_type);}
  
  //smart photon insertion
  EXTERN_PARS int use_photon_field INIT_TO(true);
  inline void read_use_photon_field()
  {read_str_int("UsePhotonField",&use_photon_field);}
  
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
  
  //compute mass or/and QED corrections
  EXTERN_PARS int compute_mass_corrections INIT_TO(true);
  EXTERN_PARS int compute_QED_corrections INIT_TO(true);
  inline void read_corrections_to_compute()
  {
    read_str_int("ComputeMassCorrections",&compute_mass_corrections);
    read_str_int("ComputeQEDCorrections",&compute_QED_corrections);
  }
  
  //number of configurations
  inline void read_ngauge_conf()
  {read_str_int("NGaugeConf",&ngauge_conf);}
  
  //number of sources
  EXTERN_PARS int nsources INIT_TO(1);
  inline void read_nsources()
  {read_str_int("NSources",&nsources);}
  
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
}

#undef EXTERN_PARS
#undef INIT_TO

#endif
