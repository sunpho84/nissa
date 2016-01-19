#ifndef _PARSER_DRIVER_HPP
#define _PARSER_DRIVER_HPP

#include "nissa.hpp"

using namespace nissa;

class driver_t
{
public:
  void *scanner;
  FILE *fin;
  driver_t(FILE *file);
  
  //geometry
  int L;
  int T;
  int def_L(){return 4;}
  int def_T(){return 4;}
  
  //walltime and seed
  int walltime;
  int seed;
  int def_walltime(){return 3500;}
  int def_seed(){return 23526211;}
  
  //gauge action
  double beta;
  gauge_action_name_t gauge_action_name;
  double def_beta(){return 6.0;}
  gauge_action_name_t def_gauge_action_name(){return WILSON_GAUGE_ACTION;}
  
  //topo potential pars
  topotential_pars_t topotential_pars;
  
  //content of quarks
  std::vector<quark_content_t> quarks;
  
  //stout parameters
  stout_pars_t stout_pars;
  
  //background em field
  em_field_pars_t em_field_pars;
  
  //fermionic measures
  nucleon_corr_meas_pars_t nucleon_corr_meas_pars;
  meson_corr_meas_pars_t meson_corr_meas_pars;
  fermionic_putpourri_meas_pars_t fermionic_putpourri_meas_pars;
  quark_rendens_meas_pars_t quark_rendens_meas_pars;
  magnetization_meas_pars_t magnetization_meas_pars;
  
  //gauge measures
  gauge_obs_meas_pars_t plaq_pol_meas_pars;
  
  //evolution and conf
  hmc_evol_pars_t evol_pars;
  conf_pars_t conf_pars;
  
  int master_fprintf_geometry(FILE *fout,int full)
  {
    int nprinted=0;
    if(full||L!=def_L()) nprinted+=nissa::master_fprintf(fout,"L\t\t=\t%d\n",L);
    if(full||T!=def_T()) nprinted+=nissa::master_fprintf(fout,"T\t\t=\t%d\n",T);
    
    return nprinted;
  }
  
  int master_fprintf_betaact(FILE *fout,int full)
  {
    int nprinted=0;
    if(full||(beta!=def_beta())) nprinted+=nissa::master_fprintf(fout,"Beta\t\t=\t%lg\n",beta);
    if(full||(gauge_action_name!=def_gauge_action_name()))
      {
	nprinted+=nissa::master_fprintf(fout,"GaugeAction\t=\t");
	switch(gauge_action_name)
	  {
	  case WILSON_GAUGE_ACTION: nprinted+=nissa::master_fprintf(fout,"Wilson");break;
	  case TLSYM_GAUGE_ACTION: nprinted+=nissa::master_fprintf(fout,"tlSym");break;
	  case IWASAKI_GAUGE_ACTION: nprinted+=nissa::master_fprintf(fout,"Iwasaki");break;
	  default:crash("unknown gauge action %d",(int)gauge_action_name);
	  }
	nprinted+=nissa::master_fprintf(fout,"\n");
      }
    
    return nprinted;
  }
  
  int master_fprintf_walltime_seed(FILE *fout,int full)
  {
    int nprinted=0;
    if(full||walltime!=def_walltime()) nprinted+=nissa::master_fprintf(fout,"Walltime\t\t=\t%d\n",walltime);
    if(full||seed!=def_seed()) nprinted+=nissa::master_fprintf(fout,"Seed\t\t\t=\t%d\n",seed);
    
    return nprinted;
  }
  
  int master_fprintf(FILE *fout,bool full=false)
  {
    int nprinted=0;
    
    //geometry
    if(master_fprintf_geometry(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    //beta and action
    if(master_fprintf_betaact(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    //topotential
    if(topotential_pars.master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    //quarks
    for(size_t i=0;i<quarks.size();i++) if(quarks[i].master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    //global stout pars
    if(full||stout_pars.nlevels!=stout_pars.def_nlevels()||stout_pars.rho!=stout_pars.def_rho()) nprinted+=nissa::master_fprintf(fout,"GlobalStoutPars\n");
    if(stout_pars.master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    //global em field pars
    if(em_field_pars.master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    //fermionic measures
    if(meson_corr_meas_pars.master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    if(nucleon_corr_meas_pars.master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    if(fermionic_putpourri_meas_pars.master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    if(quark_rendens_meas_pars.master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    if(magnetization_meas_pars.master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    //gauge masures
    if(plaq_pol_meas_pars.master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    //hmc evolution
    if(evol_pars.master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    //walltime and seed
    if(master_fprintf_walltime_seed(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    //configuration
    if(conf_pars.master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    
    return nprinted;
  }
  
private:
  driver_t() {}
  
protected:
  void init_scanner();
  void destroy_scanner();
};

int parser_parse(driver_t *driver);

#endif
