#ifndef _TOPOLOGICAL_ACTION_HPP
#define _TOPOLOGICAL_ACTION_HPP

#include <geometry/geometry_mix.hpp>
#include <io/input.hpp>
#include <measures/gauge/topological_charge.hpp>
#include <new_types/metadynamics.hpp>
#include <operations/smearing/stout.hpp>
#include <routines/ios.hpp>

namespace nissa
{
  constexpr char topo_file_name[]="topo_potential"; //improve
  
  //parameters to add topological potential
  struct topotential_pars_t :
    meta_pars_t
  {
    int flag;
    
    int def_flag() const
    {
      return 0;
    }
    
    double theta;
    
    double def_theta() const
    {
      return 0.0;
    }
    
    stout_pars_t stout_pars;
    
    void store_if_needed(const EoField<quad_su3>& ext_conf,
			 const int& iconf) const
    {
      if(flag==2 and iconf%each==0 and iconf>=after)
	{
	  CRASH("reimplement");
	  
	  // double charge;
	  // eo_ptr<quad_su3> conf;
	  // if(stout_pars.nlevels==0)
	  //   {
	  //     conf[0]=ext_conf[0];
	  //     conf[1]=ext_conf[1];
	  //   }
	  // else
	  //   {
	  //     conf[0]=nissa_malloc("stout_conf_e",locVolh+bord_volh+edge_volh,quad_su3);
	  //     conf[1]=nissa_malloc("stout_conf_o",locVolh+bord_volh+edge_volh,quad_su3);
	  //     stout_smear(conf,ext_conf,&stout_pars);
	  //   }
	  
	  // //compute topocharge
	  // total_topological_charge_eo_conf(&charge,conf);
	  // MASTER_PRINTF("Topological charge to be stored: %lg\n",charge);
	  // update(iconf,charge);
	  
	  // //free if needed
	  // if(stout_pars.nlevels!=0)
	  //   {
	  //     nissa_free(conf[0]);
	  //     nissa_free(conf[1]);
	  //   }
	}
      }
    
    int master_fprintf(FILE *fout,
		       const bool& full=false) const
    {
      return nissa::master_fprintf(fout,"%s",get_str().c_str());
    }
    
    std::string get_str(const bool& full=false) const
    {
      std::ostringstream os;
      
      const char name_known[3][10]={"NONE","ORDINARY","META"};
      
      if(full or flag!=def_flag())
	os<<"TopoPotential\t=\t"<<name_known[flag]<<"\n";
      
      switch(flag)
	{
	case 0:
	  break;
	case 1:
	  os<<" Theta\t\t=\t"<<theta<<"\n";
	  break;
	case 2:
	  os<<meta_pars_t::get_str(full);
	  os<<stout_pars.get_str(full);
	  break;
	}
      
      return os.str();
    }
    
    bool is_nonstandard() const
    {
      return
	meta_pars_t::is_nonstandard() or
	flag!=def_flag() or
	theta!=def_theta() or
	stout_pars.is_nonstandard();
    }
    
    topotential_pars_t() :
      meta_pars_t(),
      flag(def_flag()),
      theta(def_theta()){}
  };
  
  //compute the topodynamical potential using past history
  inline double topodynamical_potential(const double& Q,
					const topotential_pars_t& pars)
  {
    return pars.compute_pot(Q);
  }
  
  //draw the topodynamical potential
  inline void save_topodynamical_potential(const topotential_pars_t& pars)
  {
    pars.save(topo_file_name);
  }
  
  inline void load_topodynamical_potential(topotential_pars_t& pars,
					   const bool& mandatory)
  {
    if(file_exists(topo_file_name))
      pars.load(topo_file_name);
    else
      if(mandatory)
	CRASH("%s file not found when mandatory present",topo_file_name);
      else
	VERBOSITY_LV2_MASTER_PRINTF("%s not found, skipping reading",topo_file_name);
  }
  
  //Compute the topological action
  inline double topotential_action(const EoField<quad_su3>& conf,
				   const topotential_pars_t &pars)
  {
    //compute topocharge
    double Q;
    if(pars.stout_pars.nlevels)
      {
	EoField<quad_su3> smeConf("smeConf",WITH_HALO);
        stout_smear(smeConf,conf,pars.stout_pars);
	Q=total_topological_charge_eo_conf(smeConf);
      }
    else
      Q=total_topological_charge_eo_conf(conf);
    
    //compute according to flag
    double topo_action=0.0;
    switch(pars.flag)
      {
      case 1:
	topo_action=Q*pars.theta;
	break;
      case 2:
	topo_action=topodynamical_potential(Q,pars);
	break;
      default:
	CRASH("unknown flag %d",pars.flag);
      }
    
    return topo_action;
  }
  
  //lx version
  inline double topotential_action(const LxField<quad_su3>& lx_conf,
			    const topotential_pars_t &pars)
  {
    //allocate
    EoField<quad_su3> eo_conf("stout_conf",WITH_HALO_EDGES);
    
    //split and compute
    split_lx_vector_into_eo_parts(eo_conf,lx_conf);
    
    return topotential_action(eo_conf,pars);
  }
}


#endif
