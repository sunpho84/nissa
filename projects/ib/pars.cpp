#include <nissa.hpp>

#include "conf.hpp"

#define EXTERN_PARS
#include "pars.hpp"

namespace nissa
{
  //read common part of the input
  void read_input_preamble()
  {
    //init the grid
    read_init_grid();
    
    //Wall time
    read_str_double("WallTime",&wall_time);
    //Pure Wilson
    read_str_int("PureWilson",&pure_wilson);
    if(pure_wilson)
      {
	base=WILSON_BASE;
	nr=1;
	read_list_of_double_pairs("QKappaResidues",&nqmass,&qkappa,&residue);
      }
    else
      {
	base=MAX_TWIST_BASE;
	//Kappa
	read_str_double("Kappa",&kappa);
	//One or two r
	read_str_int("NR",&nr);
	//Masses and residue
	read_list_of_double_pairs("QMassResidues",&nqmass,&qmass,&residue);
      }
  }
  
  //read all photon pars
  void read_photon_pars()
  {
    //Zero mode subtraction
    char zero_mode_sub_str[100];
    read_str_str("ZeroModeSubtraction",zero_mode_sub_str,100);
    
    if(strncasecmp(zero_mode_sub_str,"PECIONA",100)==0) photon.zms=PECIONA;
    else
      if(strncasecmp(zero_mode_sub_str,"UNNO_ALEMANNA",100)==0) photon.zms=UNNO_ALEMANNA;
      else
	if(strncasecmp(zero_mode_sub_str,"ONLY_100",100)==0) photon.zms=ONLY_100;
	else crash("Unkwnown zero mode subtraction: %s",zero_mode_sub_str);
    
    //gauge for photon propagator
    char photon_gauge_str[100];
    read_str_str("PhotonGauge",photon_gauge_str,100);
    if(strncasecmp(photon_gauge_str,"FEYNMAN",100)==0) photon.alpha=FEYNMAN_ALPHA;
    else
      if(strncasecmp(photon_gauge_str,"LANDAU",100)==0) photon.alpha=LANDAU_ALPHA;
      else
	if(strncasecmp(photon_gauge_str,"LANDAU",100)==0) read_str_double("Alpha",&photon.alpha);
	else crash("Unkwnown photon gauge: %s",photon_gauge_str);
    
    //discretization for photon propagator
    char photon_discrete_str[100];
    read_str_str("PhotonDiscretization",photon_discrete_str,100);
    if(strncasecmp(photon_discrete_str,"WILSON",100)==0) photon.c1=WILSON_C1;
    else
      if(strncasecmp(photon_discrete_str,"TLSYM",100)==0) photon.c1=TLSYM_C1;
      else crash("Unkwnown photon discretization: %s",photon_discrete_str);
    
    //compute the tadpole summing all momentum
    compute_tadpole(tadpole,photon);
  }
}
