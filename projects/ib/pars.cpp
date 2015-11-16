#include <nissa.hpp>

#include "conf.hpp"

#define EXTERN
 #include "pars.hpp"
#undef EXTERN

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
}
