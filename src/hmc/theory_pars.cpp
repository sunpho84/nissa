#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "theory_pars.hpp"

namespace nissa
{
  int theory_pars_t::master_fprintf(FILE *fout,int full)
  {
    int nprinted=0;
    //header
    if(full||is_nonstandard()) nprinted+=nissa::master_fprintf(fout,"\nTheory\n");
    //gauge action
    if(full||(gauge_action_name!=def_gauge_action_name()))
      {
	nprinted+=nissa::master_fprintf(fout," GaugeAction\t=\t");
	switch(gauge_action_name)
	  {
	  case WILSON_GAUGE_ACTION: nprinted+=nissa::master_fprintf(fout,"Wilson");break;
	  case TLSYM_GAUGE_ACTION: nprinted+=nissa::master_fprintf(fout,"tlSym");break;
	  case IWASAKI_GAUGE_ACTION: nprinted+=nissa::master_fprintf(fout,"Iwasaki");break;
	  default:crash("unknown gauge action %d",(int)gauge_action_name);
	  }
	nprinted+=nissa::master_fprintf(fout,"\n");
      }
    //beta
    if(full||(beta!=def_beta())) nprinted+=nissa::master_fprintf(fout," Beta\t\t=\t%lg\n",beta);
    //topotential_pars
    if(topotential_pars.master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    //quarks
    for(size_t i=0;i<quarks.size();i++) if(quarks[i].master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    //stout pars
    if(stout_pars.master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    //global em field pars
    if(em_field_pars.master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    if(full||is_nonstandard()) nprinted+=nissa::master_fprintf(fout,"\n");
    
    return nprinted;
  }
}
