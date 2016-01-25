#ifndef _MULTIPSEUDO_RHMC_STEP_HPP
#define _MULTIPSEUDO_RHMC_STEP_HPP

#include "hmc/theory_pars.hpp"
#include "new_types/rat_approx.hpp"
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
    double def_pf_action_residue(){return 1e-12;}
    double def_md_residue(){return 1e-6;}
    int def_nmd_steps(){return 13;}
    int def_ngauge_substeps(){return 5;}
    
    std::vector<int> npseudo_fs;
    rat_approx_t *rat_appr;
    
    int master_fprintf(FILE *fout,bool full=false)
    {
      int nprinted=0;
      
      if(full||is_nonstandard())
	{
	  nissa::master_fprintf(fout,"Evolution\n");
	  if(full||id_sea_theory!=def_id_sea_theory()) nprinted+=nissa::master_fprintf(fout," IdSeaTheory\t=\t%d\n",id_sea_theory);
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
    
    std::string def_path(){return "conf";}
    std::string def_store_path(){return "stored_conf.08%d";}
    int def_store_each(){return 10;}
    int def_store_running(){return 1;}
    start_conf_cond_t def_start_cond(){return COLD_START_COND;}
    
    int master_fprintf(FILE *fout,bool full)
    {
      int nprinted=0;
      
      if(full||is_nonstandard())
	{
	  nissa::master_fprintf(fout,"GaugeConf\n");
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
    
    int is_nonstandard()
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
    int is_stag,npf;
    color **stag;
    spincolor **Wils;
    
    void create(int npf,int is_stag);
    void destroy();
  };
  
  double multipseudo_rhmc_step(quad_su3 **out_conf,quad_su3 **in_conf,theory_pars_t &theory_pars,
				 hmc_evol_pars_t &simul_pars,int itraj);
}

#endif
