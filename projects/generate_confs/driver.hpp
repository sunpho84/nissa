#ifndef _DRIVER_HPP
#define _DRIVER_HPP

#include "nissa.hpp"

namespace nissa
{
  //return 1 if we need to measure
  template <class T> int measure_is_due(T &pars,int iconf)
  {return (pars.each>0)&&(iconf%pars.each==0)&&(iconf>=pars.after);}
  
  class driver_t
  {
  public:
    void *scanner;
    FILE *fin;
    driver_t(FILE *file);
    ~driver_t() {destroy_scanner();}
    
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
    
    //parameters of actions
    std::vector<theory_pars_t> theories;
    int ntheories(){return theories.size();}
    theory_pars_t &sea_theory(){return theories[evol_pars.id_sea_theory];}
    
    //fermionic measures
    std::vector<meson_corr_meas_pars_t> meson_corr_meas;
    std::vector<nucleon_corr_meas_pars_t> nucleon_corr_meas;
    std::vector<fermionic_putpourri_meas_pars_t> fermionic_putpourri_meas;
    std::vector<quark_rendens_meas_pars_t> quark_rendens_meas;
    std::vector<magnetization_meas_pars_t> magnetization_meas;
    
    //check if any measure is due
    template <class T> int measure_is_due(std::vector<T> &pars,int itheory,int iconf)
    {
      int due=0;
      for(typename std::vector<T>::iterator it=pars.begin();it!=pars.end();it++) due|=it->measure_is_due(itheory,iconf);
      return due;
    }
    int any_fermionic_measure_is_due(int itheory,int iconf)
    {
      return
	measure_is_due(fermionic_putpourri_meas,itheory,iconf)||
	measure_is_due(magnetization_meas,itheory,iconf)||
	measure_is_due(nucleon_corr_meas,itheory,iconf)||
	measure_is_due(meson_corr_meas,itheory,iconf)||
	measure_is_due(quark_rendens_meas,itheory,iconf);
    }
    //print a message if a measure is due
    template <class T> bool if_meas_is_due_print(T &obj,int itheory,int iconf,const char *text)
    {
      bool a=nissa::measure_is_due(obj,iconf);
      if(a) verbosity_lv1_master_printf("Measuring %s for theory %d/%d\n",text,itheory+1,ntheories());
      return a;
    }
    
#define RANGE_FERMIONIC_MEAS(DRV,OBS)					\
    for(size_t imeas=0;imeas<DRV->NAME2(OBS,meas).size();imeas++)	\
      if(DRV->if_meas_is_due_print(DRV->NAME2(OBS,meas)[imeas],itheory,iconf,#OBS)) \
	NAME2(measure,OBS)(sme_conf,DRV->theories[itheory],DRV->NAME2(OBS,meas)[imeas],iconf,conf_created);
    
    //add
    void add_meson_corr_meas(meson_corr_meas_pars_t &m){meson_corr_meas.push_back(m);meson_corr_meas.back().itheory=ntheories()-1;}
    void add_nucleon_corr_meas(nucleon_corr_meas_pars_t &m){nucleon_corr_meas.push_back(m);nucleon_corr_meas.back().itheory=ntheories()-1;}
    void add_fermionic_putpourri_meas(fermionic_putpourri_meas_pars_t &m){fermionic_putpourri_meas.push_back(m);fermionic_putpourri_meas.back().itheory=ntheories()-1;}
    void add_quark_rendens_meas(quark_rendens_meas_pars_t &m){quark_rendens_meas.push_back(m);quark_rendens_meas.back().itheory=ntheories()-1;}
    void add_magnetization_meas(magnetization_meas_pars_t &m){magnetization_meas.push_back(m);magnetization_meas.back().itheory=ntheories()-1;}
    
    //gauge measures
    std::vector<gauge_obs_meas_pars_t> plaq_pol_meas;
    std::vector<top_meas_pars_t> top_meas;
    std::vector<poly_corr_meas_pars_t> luppoli_meas;
    std::vector<watusso_meas_pars_t> watusso_meas;
    std::vector<all_rects_meas_pars_t> all_rects_meas;
    
    //add
    void add_plaq_pol_meas(gauge_obs_meas_pars_t &m){plaq_pol_meas.push_back(m);}
    void add_top_meas(top_meas_pars_t &m){top_meas.push_back(m);}
    void add_luppoli_meas(poly_corr_meas_pars_t &m){luppoli_meas.push_back(m);}
    void add_watusso_meas(watusso_meas_pars_t &m){watusso_meas.push_back(m);watusso_meas.back();}
    void add_all_rects_meas(all_rects_meas_pars_t &m){all_rects_meas.push_back(m);all_rects_meas.back();}
    
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
      //theories
      for(size_t i=0;i<theories.size();i++) if(theories[i].master_fprintf(fout,full)) {nprinted++;nissa::master_fprintf(fout,"\n");}
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
}

#endif

