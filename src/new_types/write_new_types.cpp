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
  
  //convert a gauge action name into a str
  std::string get_gauge_action_name_str(gauge_action_name_t name)
  {
    switch(name)
      {
      case WILSON_GAUGE_ACTION:return "Wilson";break;
      case TLSYM_GAUGE_ACTION:return "tlSym";break;
      case IWASAKI_GAUGE_ACTION:return "Iwasaki";break;
      default:return "Unknown";
      }
  }
  
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
  
  //print stout_pars
  int stout_pars_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    if(full||is_nonstandard())
      {
	 nprinted+=nissa::master_fprintf(fout,"StoutPars\n");
	 if(full||nlevels!=def_nlevels()) nprinted+=nissa::master_fprintf(fout," NLevels\t=\t%d\n",nlevels);
	 if(full||rho!=def_rho()) nprinted+=nissa::master_fprintf(fout," Rho\t\t=\t%lg\n",rho);
      }
    
    return nprinted;
  }
  
  //print topotential_pars_t
  int topotential_pars_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    const char name_known[3][10]={"None","","Meta"};
    if(full||flag!=def_flag()) nprinted+=nissa::master_fprintf(fout,"TopoPotential\t=\t%s\n",name_known[flag]);
    switch(flag)
      {
      case 0:break;
      case 1:nprinted+=nissa::master_fprintf(fout,"Theta\t\t%lg\n",theta);break;
      case 2:
	nprinted+=meta_pars_t::master_fprintf(fout);
	stout_pars.master_fprintf(fout);
	break;
      }
    
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
  
  //print hmc_evol_pars_t
  int hmc_evol_pars_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    
    if(full||is_nonstandard())
      {
	nissa::master_fprintf(fout,"Evolution\n");
	if(full||ntraj_tot!=def_ntraj_tot()) nprinted+=nissa::master_fprintf(fout," NTrajTot\t=\t%d\n",ntraj_tot);
	if(full||skip_mtest_ntraj!=def_skip_mtest_ntraj()) nprinted+=nissa::master_fprintf(fout," SkipMetro\t=\t%d\n",skip_mtest_ntraj);
	if(full||traj_length!=def_traj_length()) nprinted+=nissa::master_fprintf(fout," TrajLength\t=\t%lg\n",traj_length);
	if(full||pf_action_residue!=def_pf_action_residue()) nprinted+=nissa::master_fprintf(fout," ActResidue\t=\t%lg\n",pf_action_residue);
	if(full||md_residue!=def_md_residue()) nprinted+=nissa::master_fprintf(fout," MdResidue\t=\t%lg\n",md_residue);
	if(full||nmd_steps!=def_nmd_steps()) nprinted+=nissa::master_fprintf(fout," NSteps\t\t=\t%d\n",nmd_steps);
	if(full||ngauge_substeps!=def_ngauge_substeps()) nprinted+=nissa::master_fprintf(fout," NSubSteps\t=\t%d\n",ngauge_substeps);
	if(full||npseudo_fs.size())
	  {
	    nprinted+=nissa::master_fprintf(fout," NPseudoFerms\t=\t{%d",npseudo_fs[0]);
	    for(size_t i=1;i<npseudo_fs.size();i++) nprinted+=nissa::master_fprintf(fout,",%d",npseudo_fs[i]);
	    nprinted+=nissa::master_fprintf(fout,"}\n");
	  }
      }
    
    return nprinted;
  }
  
  //parameters of configuration
  int conf_pars_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    
    if(full||is_nonstandard())
      {
	nissa::master_fprintf(fout,"Configuration\n");
	if(full||path!=def_path()) nprinted+=nissa::master_fprintf(fout," Path\t\t=\t\"%s\"\n",path.c_str());
	if(full||store_path!=def_store_path()) nprinted+=nissa::master_fprintf(fout," StorePath\t=\t\"%s\"\n",store_path.c_str());
	if(full||store_each!=def_store_each()) nprinted+=nissa::master_fprintf(fout," StoreEach\t=\t%d\n",store_each);
	if(full||store_running!=def_store_running()) nprinted+=nissa::master_fprintf(fout," StoreRunning\t=\t%d\n",store_running);
	if(full||start_cond!=def_start_cond())
	  {
	    nprinted+=nissa::master_fprintf(fout," StartCond\t=\t");
	    if(start_cond==HOT_START_COND) nprinted+=nissa::master_fprintf(fout,"HOT");
	    if(start_cond==COLD_START_COND) nprinted+=nissa::master_fprintf(fout,"COLD");
	    if(start_cond==UNSPEC_START_COND) crash("unspecified start cond");
	    nprinted+=nissa::master_fprintf(fout,"\n");
	  }
      }
    
    return nprinted;
  }
  
  //plaquette/polyakov
  int gauge_obs_meas_pars_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    
    nprinted+=nissa::master_fprintf(fout,"MeasPlaqPol\n");
    if(full||is_nonstandard())
      {
	if(each!=def_each()||full) nprinted+=nissa::master_fprintf(fout," Each\t\t=\t%d\n",each);
	if(after!=def_after()||full) nprinted+=nissa::master_fprintf(fout," After\t\t=\t%d\n",after);
	if(path!=def_path()||full) nprinted+=nissa::master_fprintf(fout," Path\t\t=\t\"%s\"\n",path.c_str());
      }
    
    return nprinted;
  }
  
  //print the smoothing parameter
  int smooth_pars_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    
    if(full||is_nonstandard())
      {
	if(full||method!=def_method()) nprinted+=nissa::master_fprintf(fout," SmoothMethod\t=\t%s\n",get_method_name().c_str());
	nprinted+=cool_pars.master_fprintf(fout,full);
	nprinted+=stout_pars.master_fprintf(fout,full);
	nprinted+=Wflow_pars.master_fprintf(fout,full);
      }
    
    return nprinted;
  }
  
  //print how to cool
  int cool_pars_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    
    if(full||is_nonstandard())
      {
	nprinted+=nissa::master_fprintf(fout,"CoolPars\n");
	if(full||gauge_action!=def_gauge_action()) nprinted+=nissa::master_fprintf(fout," Action\t\t=\t%s\n",get_gauge_action_name_str(gauge_action).c_str());
	if(full||nsteps!=def_nsteps()) nprinted+=nissa::master_fprintf(fout," NSteps\t\t=\t%d\n",nsteps);
      }
    
    return nprinted;
  }
  
  //print how to Wilson flow
  int Wflow_pars_t::master_fprintf(FILE *fout,bool full)
  {
    int nprinted=0;
    
    if(full||is_nonstandard())
      {
	nprinted+=nissa::master_fprintf(fout,"WFlowPars\n");
	if(full||T!=def_T()) nprinted+=nissa::master_fprintf(fout," FlowTime\t=\t%lg\n",T);
	if(full||dt!=def_dt()) nprinted+=nissa::master_fprintf(fout," InteStep\t=\t%lg\n",dt);
      }
    
    return nprinted;
  }
}
