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
  int quark_content_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    if(full||name!=def_name()) nprinted+=nissa::master_fprintf(fout,"Quark\t\t\"%s\"\n",name.c_str());
    if(full||deg!=def_deg()) nprinted+=nissa::master_fprintf(fout,"Degeneracy\t=\t%d\n",deg);
    if(full||mass!=def_mass()) nprinted+=nissa::master_fprintf(fout,"Mass\t\t=\t%lg\n",mass);
    if(full||re_pot!=def_re_pot()) nprinted+=nissa::master_fprintf(fout,"RePot\t\t=\t%lg\n",re_pot);
    if(full||im_pot!=def_im_pot()) nprinted+=nissa::master_fprintf(fout,"ImPot\t\t=\t%lg\n",re_pot);
    if(full||charge!=def_charge()) nprinted+=nissa::master_fprintf(fout,"ElecCharge\t=\t%lg\n",charge);

    return nprinted;
  }
  
  //print stout_pars
  int stout_pars_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    if(full||nlevels!=def_nlevels()) nprinted+=nissa::master_fprintf(fout,"NLevels\t\t=\t%d\n",nlevels);
    if(full||rho!=def_rho()) nprinted+=nissa::master_fprintf(fout,"Rho\t\t=\t%lg\n",rho);

    return nprinted;
  }
  
  //print topotential_pars_t
  int topotential_pars_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
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
    
    return nprinted;
  }
  
  //print em_field_pars
  int em_field_pars_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    nissa::master_fprintf(fout,"BkgrdEMField");
    if(flag)
      {
	nissa::master_fprintf(fout,"\n");
	if(full||fabs(E[0])>1e-14) nprinted+=nissa::master_fprintf(fout,"Ex\t\t=\t%lg\n",E[0]);
	if(full||fabs(E[1])>1e-14) nprinted+=nissa::master_fprintf(fout,"Ey\t\t=\t%lg\n",E[1]);
	if(full||fabs(E[2])>1e-14) nprinted+=nissa::master_fprintf(fout,"Ez\t\t=\t%lg\n",E[2]);
	if(full||fabs(B[0])>1e-14) nprinted+=nissa::master_fprintf(fout,"Bx\t\t=\t%lg\n",B[0]);
	if(full||fabs(B[1])>1e-14) nprinted+=nissa::master_fprintf(fout,"By\t\t=\t%lg\n",B[1]);
	if(full||fabs(B[2])>1e-14) nprinted+=nissa::master_fprintf(fout,"Bz\t\t=\t%lg\n",B[2]);
      }
    else nprinted+=nissa::master_fprintf(fout,"\t\tNo\n");
    
    return nprinted;
  }
  
  //print hmc_evol_pars_t
  int hmc_evol_pars_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    
    return nprinted;
  }
  
  //pseudo correlators
  int pseudo_corr_meas_pars_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    
    if(flag||full)
      {
	nprinted+=nissa::master_fprintf(fout,"PseudoCorrelators\n");
	if(flag!=1||full) nprinted+=nissa::master_fprintf(fout,"Each\t\t=\t%d\n",flag);
	if(path!=def_path()||full) nprinted+=nissa::master_fprintf(fout,"Path\t\t=\t\"%s\"\n",path.c_str());
	if(residue!=def_residue()||full) nprinted+=nissa::master_fprintf(fout,"Residue\t\t=\t\"%lg\"\n",residue);
	if(nhits!=def_nhits()||full) nprinted+=nissa::master_fprintf(fout,"NHits\t\t=\t%d\n",nhits);
      }
    else if(full) nprinted+=nissa::master_fprintf(fout,"PseudoCorrelators No\n");
    
    return nprinted;
  }
  
  //pseudo putpourri
  int fermionic_putpourri_meas_pars_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    
    if(flag||full)
      {
	nprinted+=nissa::master_fprintf(fout,"FermionicPutpourri\n");
	if(flag!=1||full) nprinted+=nissa::master_fprintf(fout,"Each\t\t=\t%d\n",flag);
	if(path!=def_path()||full) nprinted+=nissa::master_fprintf(fout,"Path\t\t=\t\"%s\"\n",path.c_str());
	if(residue!=def_residue()||full) nprinted+=nissa::master_fprintf(fout,"Residue\t\t=\t\"%lg\"\n",residue);
	if(compute_susceptivities!=def_compute_susceptivities()||full) nprinted+=nissa::master_fprintf(fout,"Residue\t\t=\t\"%lg\"\n",compute_susceptivities);
	if(ncopies!=def_ncopies()||full) nprinted+=nissa::master_fprintf(fout,"NCopies\t\t=\t%d\n",ncopies);
	if(nhits!=def_nhits()||full) nprinted+=nissa::master_fprintf(fout,"NHits\t\t=\t%d\n",nhits);
      }
    else if(full) nprinted+=nissa::master_fprintf(fout,"FermionicPutpourri No\n");
    
    return nprinted;
  }
  
  //rendens
  int quark_rendens_meas_pars_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    
    if(flag||full)
      {
	nprinted+=nissa::master_fprintf(fout,"QuarkRendens\n");
	if(flag!=1||full) nprinted+=nissa::master_fprintf(fout,"Each\t\t=\t%d\n",flag);
	if(path!=def_path()||full) nprinted+=nissa::master_fprintf(fout,"Path\t\t=\t\"%s\"\n",path.c_str());
	if(residue!=def_residue()||full) nprinted+=nissa::master_fprintf(fout,"Residue\t\t=\t\"%lg\"\n",residue);
	if(nhits!=def_nhits()||full) nprinted+=nissa::master_fprintf(fout,"NHits\t\t=\t%d\n",nhits);
      }
    else if(full) nprinted+=nissa::master_fprintf(fout,"QuarkRendens No\n");
    
    return nprinted;
  }
}
