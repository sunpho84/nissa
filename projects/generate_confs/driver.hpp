#ifndef _DRIVER_HPP
#define _DRIVER_HPP

#include "nissa.hpp"

using namespace nissa;

class driver_t
{
public:
  void *scanner;
  FILE *fin;
  driver_t(FILE *file);
  
  //geometry
  int LX,LY,LZ;
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
  std::vector<meson_corr_meas_pars_t> meson_corr_meas;
  std::vector<nucleon_corr_meas_pars_t> nucleon_corr_meas;
  std::vector<fermionic_putpourri_meas_pars_t> fermionic_putpourri_meas;
  std::vector<quark_rendens_meas_pars_t> quark_rendens_meas;
  std::vector<magnetization_meas_pars_t> magnetization_meas;
  
  //gauge measures
  std::vector<gauge_obs_meas_pars_t> plaq_pol_meas;
  std::vector<top_meas_pars_t> top_meas;
  std::vector<poly_corr_meas_pars_t> luppoli_meas;
  std::vector<watusso_meas_pars_t> watusso_meas;
  std::vector<all_rect_meas_pars_t> all_rect_meas;
  
  //mode of running
  enum run_mode_t{EVOLUTION_MODE,ANALYSIS_MODE};
  run_mode_t run_mode;
  run_mode_t def_run_mode(){return EVOLUTION_MODE;}
  
  //evolution and conf
  hmc_evol_pars_t evol_pars;
  conf_pars_t conf_pars;
  std::vector<std::string> an_conf_list;
  
  int master_fprintf_geometry(FILE *fout,int full)
  {
    int nprinted=0;
    int fX=(LX!=def_L()),fY=(LY!=def_L()),fZ=(LZ!=def_L());
    int fL=(LX==LY&&LX==LZ),fT=(T!=def_T());;
    if(full||fX||fY||fZ||fT) nprinted+=nissa::master_fprintf(fout,"Geometry\n");
    if(full||fT) nprinted+=nissa::master_fprintf(fout," T\t\t=\t%d\n",T);
    if(!fL)
      {
	if(full||fX) nprinted+=nissa::master_fprintf(fout," LX\t\t=\t%d\n",LX);
	if(full||fY) nprinted+=nissa::master_fprintf(fout," LY\t\t=\t%d\n",LY);
	if(full||fZ) nprinted+=nissa::master_fprintf(fout," LZ\t\t=\t%d\n",LZ);
      }
    else if(full||(fX||fY||fZ)) nprinted+=nissa::master_fprintf(fout," L\t\t=\t%d\n",LX);
    
    return nprinted;
  }
  
  int master_fprintf_betaact(FILE *fout,int full)
  {
    int nprinted=0;
    if(full||(beta!=def_beta())||(gauge_action_name!=def_gauge_action_name())) nprinted+=nissa::master_fprintf(fout,"GaugePars\n");
    if(full||(beta!=def_beta())) nprinted+=nissa::master_fprintf(fout," Beta\t\t=\t%lg\n",beta);
    if(full||(gauge_action_name!=def_gauge_action_name()))
      {
	nprinted+=nissa::master_fprintf(fout," Action\t\t=\t");
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
    if(full||walltime!=def_walltime()||seed!=def_seed()) nprinted+=nissa::master_fprintf(fout,"Run\n");
    if(full||walltime!=def_walltime()) nprinted+=nissa::master_fprintf(fout," Walltime\t=\t%d\n",walltime);
    if(full||seed!=def_seed()) nprinted+=nissa::master_fprintf(fout," Seed\t\t=\t%d\n",seed);
    
    return nprinted;
  }
  
  //print a whole vector
  template <class T> int master_printf_vector(FILE *fout,std::vector<T> &v,int full)
  {
    int nprinted=0;
    for(typename std::vector<T>::iterator it=v.begin();it!=v.end();it++) if(it->master_fprintf(fout,full)){nprinted+=nissa::master_fprintf(fout,"\n");}
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
    if(stout_pars.master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    //global em field pars
    if(em_field_pars.master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    //fermionic measures
    nprinted+=master_printf_vector(fout,meson_corr_meas,full);
    nprinted+=master_printf_vector(fout,nucleon_corr_meas,full);
    nprinted+=master_printf_vector(fout,fermionic_putpourri_meas,full);
    nprinted+=master_printf_vector(fout,quark_rendens_meas,full);
    nprinted+=master_printf_vector(fout,magnetization_meas,full);
    //gauge masures
    nprinted+=master_printf_vector(fout,plaq_pol_meas,full);
    nprinted+=master_printf_vector(fout,top_meas,full);
    //walltime and seed
    if(master_fprintf_walltime_seed(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
    
    switch(run_mode)
      {
      case EVOLUTION_MODE:
	if(evol_pars.master_fprintf(fout,full)) nprinted+=nissa::master_fprintf(fout,"\n");
	if(conf_pars.master_fprintf(fout,full)) nprinted+=nissa::master_fprintf(fout,"\n");
	break;
      case ANALYSIS_MODE:
	if(full||def_run_mode()!=ANALYSIS_MODE) nprinted+=nissa::master_fprintf(fout,"Analysis\n");
	if(an_conf_list.size())
	  {
	    nprinted+=nissa::master_fprintf(fout," ConfList\t=\t{\"%s\"",an_conf_list[0].c_str());
	    for(size_t i=1;i<an_conf_list.size();i++) nprinted+=nissa::master_fprintf(fout,",\"%s\"",an_conf_list[i].c_str());
	    nprinted+=nissa::master_fprintf(fout,"}\n");
	  }
	break;
      }
    
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
