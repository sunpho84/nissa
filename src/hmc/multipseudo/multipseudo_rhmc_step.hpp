#ifndef _MULTIPSEUDO_RHMC_STEP_HPP
#define _MULTIPSEUDO_RHMC_STEP_HPP

#include "base/random.hpp"
#include "hmc/theory_pars.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/rat_approx.hpp"
#include "geometry/geometry_eo.hpp"
#include "operations/gaugeconf.hpp"

namespace nissa
{
  //evolution parameters for hybrid monte carlo
  struct hmc_evol_pars_t
  {
    int id_sea_theory;
    int ntraj_tot;
    int skip_mtest_ntraj;
    double traj_length;
    double pf_action_residue;
    double md_residue;
    int nmd_steps;
    int ngauge_substeps;
    
    int def_id_sea_theory(){return 0;}
    int def_ntraj_tot(){return 100;}
    int def_skip_mtest_ntraj(){return 30;}
    double def_traj_length(){return 1.0;}
    double def_pf_action_residue(){return 1e-16;}
    double def_md_residue(){return 1e-8;}
    int def_nmd_steps(){return 11;}
    int def_ngauge_substeps(){return 5;}
    
    std::vector<int> npseudo_fs;
    
    int master_fprintf(FILE *fout,int full) {return nissa::master_fprintf(fout,"%s",get_str().c_str());}
    std::string get_str(int full=false)
    {
      std::ostringstream os;
      
      if(full||is_nonstandard())
	{
	  os<<"Evolution\n";
	  if(full||id_sea_theory!=def_id_sea_theory()) os<<" IdSeaTheory\t=\t"<<id_sea_theory<<"\n";
	  if(full||ntraj_tot!=def_ntraj_tot()) os<<" NTrajTot\t=\t"<<ntraj_tot<<"\n";
	  if(full||skip_mtest_ntraj!=def_skip_mtest_ntraj()) os<<" SkipMetro\t=\t"<<skip_mtest_ntraj<<"\n";
	  if(full||traj_length!=def_traj_length()) os<<" TrajLength\t=\t"<<traj_length<<"\n";
	  if(full||pf_action_residue!=def_pf_action_residue()) os<<" ActResidue\t=\t"<<pf_action_residue<<"\n";
	  if(full||md_residue!=def_md_residue()) os<<" MdResidue\t=\t"<<md_residue<<"\n";
	  if(full||nmd_steps!=def_nmd_steps()) os<<" NSteps\t\t=\t"<<nmd_steps<<"\n";
	  if(full||ngauge_substeps!=def_ngauge_substeps()) os<<" NSubSteps\t=\t"<<ngauge_substeps<<"\n";
	  if(full||npseudo_fs.size())
	    {
	      os<<" NPseudoFerms\t=\t{";
	      if(npseudo_fs.size()) os<<npseudo_fs[0];
	      for(size_t i=1;i<npseudo_fs.size();i++) os<<","<<npseudo_fs[i];
	      os<<"}\n";
	    }
	}
      
      return os.str();
    }
    
    int is_nonstandard()
    {
      return
	id_sea_theory!=def_id_sea_theory()||
	ntraj_tot!=def_ntraj_tot()||
	skip_mtest_ntraj!=def_skip_mtest_ntraj()||
	traj_length!=def_traj_length()||
	pf_action_residue!=def_pf_action_residue()||
	md_residue!=def_md_residue()||
	nmd_steps!=def_nmd_steps()||
	ngauge_substeps!=def_ngauge_substeps();
    }
    
    hmc_evol_pars_t() :
      id_sea_theory(def_id_sea_theory()),
      ntraj_tot(def_ntraj_tot()),
      skip_mtest_ntraj(def_skip_mtest_ntraj()),
      traj_length(def_traj_length()),
      pf_action_residue(def_pf_action_residue()),
      md_residue(def_md_residue()),
      nmd_steps(def_nmd_steps()),
      ngauge_substeps(def_ngauge_substeps()) {}
  };
  
  struct conf_pars_t
  {
    std::string path;
    std::string store_path;
    int store_each;
    int store_running;
    start_conf_cond_t start_cond;
    
    std::string def_path() const
    {
      return "conf";
    }
    
    std::string def_store_path() const
    {
      return "stored_conf.%08d";
    }
    
    int def_store_each() const
    {
      return 10;
    }
    
    int def_store_running() const
    {
      return 1;
    }
    
    start_conf_cond_t def_start_cond() const
    {
      return COLD_START_COND;
    }
        
    int master_fprintf(FILE *fout,
		       const bool& full) const
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool full=false) const
    {
      std::ostringstream os;
      
      if(full||is_nonstandard())
	{
	  os<<"GaugeConf\n";
	  if(full||path!=def_path()) os<<" Path\t\t=\t\""<<path<<"\"\n";
	  if(full||store_path!=def_store_path()) os<<" StorePath\t=\t\""<<store_path<<"\"\n";
	  if(full||store_each!=def_store_each()) os<<" StoreEach\t=\t"<<store_each<<"\n";
	  if(full||store_running!=def_store_running()) os<<" StoreRunning\t=\t"<<store_running<<"\n";
	  if(full||start_cond!=def_start_cond())
	    {
	      os<<" StartCond\t=\t";
	      if(start_cond==HOT_START_COND) os<<"HOT";
	      if(start_cond==COLD_START_COND) os<<"COLD";
	      if(start_cond==UNSPEC_START_COND) crash("unspecified start cond");
	      os<<"\n";
	    }
	}
      
      return os.str();
    }
    
    int is_nonstandard() const
    {
      return
	path!=def_path()||
	store_path!=def_store_path()||
	store_each!=def_store_each()||
	store_running!=def_store_running()||
	start_cond!=def_start_cond();
    }
    
    conf_pars_t() :
      path(def_path()),
      store_path(def_store_path()),
      store_each(def_store_each()),
      store_running(def_store_running()),
      start_cond(def_start_cond()) {}
  };
  
  /////////////////////////////
  
  //to hold rhmc fields
  struct pseudofermion_t
  {
    int is_stag;
    // double *double_ptr;
    
    EvnField<color> *stag;
    EvnField<spincolor> *Wils;
    
    //fill to random
    void fill(const rnd_t& rtype=RND_GAUSS,
	      const int& twall=-1,
	      const int& par=EVN,
	      const int& dir=0)
    {
      if(is_stag) generate_fully_undiluted_eo_source(stag->castFieldCoverage<EVEN_OR_ODD_SITES,true>(),rtype,twall,par,dir);
      else generate_fully_undiluted_eo_source(Wils->castFieldCoverage<EVEN_OR_ODD_SITES,true>(),rtype,twall,par,dir);
    }
    
    // //normalize and return the size
    // double normalize(const pseudofermion_t& other_vector,
    // 		     const double& norm=1)
    // {
    //   double other_norm;
    //   double_vector_normalize(&other_norm,double_ptr,other_vector.double_ptr,norm,ndoubles);
      
    //   return other_norm;
    // }
    
    // double normalize(const double& norm=1)
    // {
    //   const double f=1/sqrt(this->norm2());
      
    //   if(is_stag) (*stag)=f;
    //   else (*Wils)*=f;
      
    //   return 
    // }
    
    //scalar product with another vector
    double scal_prod_with(const pseudofermion_t& oth) const
    {
      double res;
      
      if(is_stag)
	res=stag->realPartOfScalarProdWith(*oth.stag);
      else
	res=Wils->realPartOfScalarProdWith(*oth.Wils);
      
      return res;
    }
    
    //return the squared norm
    double norm2() const
    {
      return scal_prod_with(*this);
    }
    
    //allocate and mark size
    void create(const ferm_discretiz::name_t& discretiz,
		const char *name="pf")
    {
      is_stag=ferm_discretiz::is_stag(discretiz);
      
      if(is_stag)
	stag=new EvnField<color>("stag",WITH_HALO);
      else
	Wils=new EvnField<spincolor>("Wils",WITH_HALO);
    }

    pseudofermion_t(const ferm_discretiz::name_t regul,
		    const char *name="pf")
    {
      create(regul,name);
    }
    
    pseudofermion_t() :
      stag(nullptr),
      Wils{nullptr}
    {
    }
    
    ~pseudofermion_t()
    {
      destroy();
    }
    
    void destroy()
    {
      if(stag) delete stag;
      if(Wils) delete Wils;
    }
  };
  
  double multipseudo_rhmc_step(EoField<quad_su3>& out_conf,
			       const EoField<quad_su3>& in_conf,
			       theory_pars_t &theory_pars,
			       hmc_evol_pars_t &simul_pars,
			       std::vector<rat_approx_t> &rat_appr,
			       const int itraj);
}

#endif
