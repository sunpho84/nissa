#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/vectors.hpp"
#include "io/input.hpp"
#include "new_types_definitions.hpp"
#include "operations/stag/nucleon.hpp"

namespace nissa
{
  int master_fprintf(FILE *stream,const char *format,...);
  
  //print quark content
  void quark_content_t::master_fprintf(FILE *fout)
  {
    nissa::master_fprintf(fout,"Quark\t\t\"%s\"\n",name.c_str());
    nissa::master_fprintf(fout,"Degeneracy\t=\t%d\n",deg);
    nissa::master_fprintf(fout,"Mass\t\t=\t%lg\n",mass);
    nissa::master_fprintf(fout,"RePot\t\t=\t%lg\n",re_pot);
    nissa::master_fprintf(fout,"ImPot\t\t=\t%lg\n",re_pot);
    nissa::master_fprintf(fout,"ElecCharge\t=\t%lg\n",charge);
  }
  
  //print stout_pars
  void stout_pars_t::master_fprintf(FILE *fout)
  {
    nissa::master_fprintf(fout,"NLevels\t\t=\t%d\n",nlevels);
    nissa::master_fprintf(fout,"Rho\t\t=\t%lg\n",rho[0][0]);
  }
  
  //print topotential_pars_t
  void topotential_pars_t::master_fprintf(FILE *fout)
  {
    const char name_known[3][10]={"None","","Meta"};
    nissa::master_fprintf(fout,"TopoPotential\t\t%s\n",name_known[flag]);
    switch(flag)
      {
      case 0:break;
      case 1:nissa::master_fprintf(fout,"Theta\t\t%lg\n",theta);break;
      case 2:
	meta_pars_t::master_fprintf(fout);
	stout_pars.master_fprintf(fout);
	break;
      }
  }
  
  //print em_field_pars
  void em_field_pars_t::master_fprintf(FILE *fout)
  {
    nissa::master_fprintf(fout,"BkgrdEMField");
    if(flag)
      {
	nissa::master_fprintf(fout,"\n");
	nissa::master_fprintf(fout,"Ex\t\t=\t%lg\n",E[0]);
	nissa::master_fprintf(fout,"Ey\t\t=\t%lg\n",E[1]);
	nissa::master_fprintf(fout,"Ez\t\t=\t%lg\n",E[2]);
	nissa::master_fprintf(fout,"Bx\t\t=\t%lg\n",B[0]);
	nissa::master_fprintf(fout,"By\t\t=\t%lg\n",B[1]);
	nissa::master_fprintf(fout,"Bz\t\t=\t%lg\n",B[2]);
      }
    else nissa::master_fprintf(fout,"\t\tNo\n");
  }
  
  //print hmc_evol_pars_t
  void hmc_evol_pars_t::master_fprintf(FILE *fout)
  {
  }
}
