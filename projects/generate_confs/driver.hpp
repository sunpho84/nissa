#ifndef _DRIVER_HPP
#define _DRIVER_HPP

#include "nissa.hpp"

namespace nissa
{
  //return 1 if we need to measure
  template <class T> int measure_is_due(T &pars,int iconf)
  {iconf-=pars.after;return (pars.each>0)&&(iconf%pars.each==0)&&(iconf>=0);}
  
  class driver_t
  {
  public:
    void *scanner;
    FILE *fin;
    driver_t(FILE *file);
    ~driver_t() {destroy_scanner();}
    
    //tag
    std::string tag;
    std::string def_tag(){return "";}
    
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
    theory_pars_t &sea_theory(){return theories[hmc_evol_pars.id_sea_theory];}
    
    //fermionic measures
    std::vector<meson_corr_meas_pars_t> meson_corr_meas;
    std::vector<nucleon_corr_meas_pars_t> nucleon_corr_meas;
    std::vector<fermionic_putpourri_meas_pars_t> fermionic_putpourri_meas;
    std::vector<quark_rendens_meas_pars_t> quark_rendens_meas;
    std::vector<chir_zumba_meas_pars_t> chir_zumba_meas;
    std::vector<spinpol_meas_pars_t> spinpol_meas;
    std::vector<qed_corr_meas_pars_t> qed_corr_meas;
    std::vector<magnetization_meas_pars_t> magnetization_meas;
    std::vector<minmax_eigenvalues_meas_pars_t> minmax_eigenvalues_meas;
    std::vector<spectr_proj_meas_pars_t> spectral_proj_meas;
    std::vector<tm_tuning_meas_pars_t> tm_tuning_meas;
    std::vector<ellesettete_meas_pars_t> ellesettete_meas;
    
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
	measure_is_due(fermionic_putpourri_meas,itheory,iconf) or
	measure_is_due(magnetization_meas,itheory,iconf) or
	measure_is_due(minmax_eigenvalues_meas,itheory,iconf) or
	measure_is_due(nucleon_corr_meas,itheory,iconf) or
	measure_is_due(meson_corr_meas,itheory,iconf) or
	measure_is_due(quark_rendens_meas,itheory,iconf) or
	measure_is_due(chir_zumba_meas,itheory,iconf) or
	measure_is_due(spinpol_meas,itheory,iconf) or
	measure_is_due(qed_corr_meas,itheory,iconf) or
	measure_is_due(spectral_proj_meas,itheory,iconf) or
	measure_is_due(tm_tuning_meas,itheory,iconf);
    }
    //print a message if a measure is due
    template <class T> bool if_meas_is_due_print(T &obj,int itheory,int iconf,const char *text)
    {
      bool a=nissa::measure_is_due(obj,iconf);
      bool b=obj.itheory==itheory;
      if(a and b) verbosity_lv1_master_printf("Measuring %s for theory %d/%d\n",text,itheory+1,ntheories());
      return a and b;
    }
    
    //add
    void add_meson_corr_meas(meson_corr_meas_pars_t &m){meson_corr_meas.push_back(m);}
    void add_nucleon_corr_meas(nucleon_corr_meas_pars_t &m){nucleon_corr_meas.push_back(m);}
    void add_fermionic_putpourri_meas(fermionic_putpourri_meas_pars_t &m){fermionic_putpourri_meas.push_back(m);}
    void add_quark_rendens_meas(quark_rendens_meas_pars_t &m){quark_rendens_meas.push_back(m);}
    void add_chir_zumba_meas(chir_zumba_meas_pars_t &m){chir_zumba_meas.push_back(m);}
    void add_spinpol_meas(spinpol_meas_pars_t &m){spinpol_meas.push_back(m);}
    void add_qed_corr_meas(qed_corr_meas_pars_t &m){qed_corr_meas.push_back(m);}
    void add_magnetization_meas(magnetization_meas_pars_t &m){magnetization_meas.push_back(m);}
    void add_minmax_eigenvalues_meas(minmax_eigenvalues_meas_pars_t &m){minmax_eigenvalues_meas.push_back(m);}
    void add_spectr_proj_meas(spectr_proj_meas_pars_t &m){spectral_proj_meas.push_back(m);}
    void add_tm_tuning_meas(tm_tuning_meas_pars_t &m){tm_tuning_meas.push_back(m);}
    void add_ellesettete_meas(ellesettete_meas_pars_t &m){ellesettete_meas.push_back(m);}

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
    bool force_unquenched;
    bool def_force_unquenched(){return false;}
    
    hmc_evol_pars_t hmc_evol_pars;
    pure_gauge_evol_pars_t quenched_evol_pars;
    conf_pars_t conf_pars;
    std::vector<std::string> an_conf_list;
    
    int master_fprintf_geometry(FILE *fout,bool full){return nissa::master_fprintf(fout,"%s",geo_get_str(full).c_str());}
    std::string geo_get_str(bool full)
    {
      std::ostringstream os;
      int fX=(LX!=def_L()),fY=(LY!=def_L()),fZ=(LZ!=def_L());
      int fL=(LX==LY&&LX==LZ),fT=(T!=def_T());;
      if(full||fX||fY||fZ||fT) os<<"Geometry\n";
      if(full||fT) os<<" T\t\t=\t"<<T<<"\n";
      if(!fL)
	{
	  if(full||fX) os<<" LX\t\t=\t"<<LX<<"\n";
	  if(full||fY) os<<" LY\t\t=\t"<<LY<<"\n";
	  if(full||fZ) os<<" LZ\t\t=\t"<<LZ<<"\n";
	}
      else if(full||(fX||fY||fZ)) os<<" L\t\t=\t"<<LX<<"\n";
      
      return os.str();
    }
    
    int master_fprintf_walltime_seed(FILE *fout,bool full){return nissa::master_fprintf(fout,"%s",walltime_seed_get_str(full).c_str());}
    std::string walltime_seed_get_str(bool full)
    {
      std::ostringstream os;
      if(full||walltime!=def_walltime()||seed!=def_seed()) os<<"Run\n";
      if(full||walltime!=def_walltime()) os<<" Walltime\t=\t"<<walltime<<"\n";
      if(full||seed!=def_seed()) os<<" Seed\t\t=\t"<<seed<<"\n";
      
      return os.str();
    }
    
    //print a whole vector
    template <class T> std::string vector_get_str(std::vector<T> &v,bool full)
    {
      std::ostringstream os;
      for(typename std::vector<T>::iterator it=v.begin();it!=v.end();it++)
	os<<it->get_str(full)<<"\n";
      return os.str();
    }
    
    //print a whole vector
    template <class T> int master_fprintf_vector(FILE *fout,std::vector<T> &v,bool full)
    {return nissa::master_fprintf(fout,"%s",get_str(full).c_str());}
    
    int master_fprintf(FILE *fout,bool full=false) {return nissa::master_fprintf(fout,"%s",get_str(full).c_str());}
    std::string get_str(bool full=false)
    {
      std::ostringstream os;
      
      //tag
      if(tag!=def_tag()||full) os<<"Tag\t\t=\t\""<<tag.c_str()<<"\"\n\n";
      //geometry
      os<<geo_get_str(full)<<"\n";
      //theories
      for(size_t i=0;i<theories.size();i++) os<<theories[i].get_str(full)<<"\n";
      //fermionic measures
      os<<vector_get_str(meson_corr_meas,full);
      os<<vector_get_str(nucleon_corr_meas,full);
      os<<vector_get_str(fermionic_putpourri_meas,full);
      os<<vector_get_str(quark_rendens_meas,full);
      os<<vector_get_str(chir_zumba_meas,full);
      os<<vector_get_str(spinpol_meas,full);
      os<<vector_get_str(qed_corr_meas,full);
      os<<vector_get_str(magnetization_meas,full);
      os<<vector_get_str(minmax_eigenvalues_meas,full);
      os<<vector_get_str(spectral_proj_meas,full);
      os<<vector_get_str(tm_tuning_meas,full);
      os<<vector_get_str(ellesettete_meas,full);
      //gauge masures
      os<<vector_get_str(plaq_pol_meas,full);
      os<<vector_get_str(top_meas,full);
      os<<vector_get_str(all_rects_meas,full);
      //walltime and seed
      os<<walltime_seed_get_str(full)<<"\n";
      
      switch(run_mode)
	{
	case EVOLUTION_MODE:
	  if(full||force_unquenched!=def_force_unquenched()) os<<"ForceUnquenched\t="<<force_unquenched<<"\n";
	  os<<hmc_evol_pars.get_str(full)<<"\n";
	  os<<quenched_evol_pars.get_str(full)<<"\n";
	  os<<conf_pars.get_str(full)<<"\n";
	  break;
	case ANALYSIS_MODE:
	  if(full||def_run_mode()!=ANALYSIS_MODE) os<<"Analysis\n";
	  if(an_conf_list.size())
	    {
	      os<<" ConfList\t=\t{\""<<an_conf_list[0].c_str()<<"\"";
	      for(size_t i=1;i<an_conf_list.size();i++) os<<",\""<<an_conf_list[i].c_str()<<"\"";
	      os<<"}\n";
	    }
	  break;
	}
      
      return os.str();
    }
    
  private:
    driver_t() {}
    
  protected:
    void init_scanner();
    void destroy_scanner();
  };
}

int parser_parse(nissa::driver_t *driver);

#endif

