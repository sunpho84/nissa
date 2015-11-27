#ifndef _PARS_HPP
#define _PARS_HPP

#include <nissa.hpp>

#ifndef EXTERN
 #define EXTERN extern
#endif

namespace nissa
{
  //convention on gospel
  const int follow_chris=0,follow_nazario=1;
  
  //define types of quark propagator used
  const int nins_kind=5;
  enum insertion_t{                    ORIGINAL,  SCALAR,  PSEUDO,  PHOTON,   TADPOLE};//,  VECTOR};
  const char ins_name[nins_kind][20]={"ORIGINAL","SCALAR","PSEUDO","PHOTON", "TADPOLE"};//, "VECTOR"};
  
  //sign of the lepton momentum
  const int norie=2;
  const int sign_orie[2]={-1,+1};
  
  const int nins=3;
  const int nlins=2;
  const int nrev=2;
  
  EXTERN int pure_wilson;
  EXTERN tm_basis_t base;
  EXTERN double kappa;
  
  EXTERN int nqmass,nr;
  EXTERN double *qmass,*qkappa,*residue;
  
  EXTERN gauge_info photon;
  EXTERN double tadpole[4];
  
  void read_input_preamble();
  void read_photon_pars();
  
  //initialize the random generator with the read seed
  inline void read_seed_start_random()
  {
    int seed;
    read_str_int("Seed",&seed);
    start_loc_rnd_gen(seed);
  }
  
  //flag to simulate in the free theory
  EXTERN int free_theory;
  inline void read_free_theory_flag()
  {read_str_int("FreeTheory",&free_theory);}
  
  //flag to make the muon with or without the external line
  EXTERN int follow_chris_or_nazario;
  inline void read_gospel_convention()
  {read_str_int("FollowChrisOrNazario",&follow_chris_or_nazario);}
  
  //noise type
  EXTERN int noise_type;
  inline void read_noise_type()
  {read_str_int("NoiseType",&noise_type);}
  
  //perform a random gauge transformation
  EXTERN int rnd_gauge_transform;
  inline void read_random_gauge_transform()
  {read_str_int("RandomGaugeTransform",&rnd_gauge_transform);}
  
  //local pion or muon current?
  EXTERN int loc_pion_curr;
  inline void read_loc_pion_curr()
  {read_str_int("LocPionCurr",&loc_pion_curr);}
  EXTERN int loc_muon_curr;
  inline void read_loc_muon_curr()
  {read_str_int("LocMuonCurr",&loc_muon_curr);}
  
  //number of configurations
  inline void read_ngauge_conf()
  {read_str_int("NGaugeConf",&ngauge_conf);}
  
  //number of sources
  EXTERN int nsources;
  inline void read_nsources()
  {read_str_int("NSources",&nsources);}
}

#undef EXTERN

#endif
