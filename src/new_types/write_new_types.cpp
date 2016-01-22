#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "io/input.hpp"
#include "operations/stag/nucleon.hpp"

namespace nissa
{
  int master_fprintf(FILE *stream,const char *format,...);
  
  //print quark content
  int quark_content_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    nprinted+=nissa::master_fprintf(fout,"Quark\t\t=\t\"%s\"\n",name.c_str());
    if(full||deg!=def_deg()) nprinted+=nissa::master_fprintf(fout," Degeneracy\t=\t%d\n",deg);
    if(full||mass!=def_mass()) nprinted+=nissa::master_fprintf(fout," Mass\t\t=\t%lg\n",mass);
    if(full||re_pot!=def_re_pot()) nprinted+=nissa::master_fprintf(fout," RePotCh\t=\t%lg\n",re_pot);
    if(full||im_pot!=def_im_pot()) nprinted+=nissa::master_fprintf(fout," ImPotCh\t=\t%lg\n",re_pot);
    if(full||charge!=def_charge()) nprinted+=nissa::master_fprintf(fout," ElecCharge\t=\t%lg\n",charge);
    
    return nprinted;
  }
  
  //print em_field_pars
  int em_field_pars_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    if(full||flag||is_nonstandard())
      {
	nissa::master_fprintf(fout,"BkgrdEMField\n");
	if(full||fabs(E[0])>1e-14) nprinted+=nissa::master_fprintf(fout," Ex\t\t=\t%lg\n",E[0]);
	if(full||fabs(E[1])>1e-14) nprinted+=nissa::master_fprintf(fout," Ey\t\t=\t%lg\n",E[1]);
	if(full||fabs(E[2])>1e-14) nprinted+=nissa::master_fprintf(fout," Ez\t\t=\t%lg\n",E[2]);
	if(full||fabs(B[0])>1e-14) nprinted+=nissa::master_fprintf(fout," Bx\t\t=\t%lg\n",B[0]);
	if(full||fabs(B[1])>1e-14) nprinted+=nissa::master_fprintf(fout," By\t\t=\t%lg\n",B[1]);
	if(full||fabs(B[2])>1e-14) nprinted+=nissa::master_fprintf(fout," Bz\t\t=\t%lg\n",B[2]);
      }
    
    return nprinted;
  }
}
