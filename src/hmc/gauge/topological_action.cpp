#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <algorithm>

#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_mix.hpp"
#include "io/input.hpp"
#include "new_types/su3.hpp"
#include "operations/smearing/stout.hpp"
#include "measures/gauge/topological_charge.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  const char topo_file_name[]="topo_potential";
  
  //compute the topodynamical potential using past history
  double topodynamical_potential(const double& Q,
				 const topotential_pars_t& pars)
  {
    return pars.compute_pot(Q);
  }

  //draw the topodynamical potential
  void save_topodynamical_potential(const topotential_pars_t& pars)
  {
    pars.save(topo_file_name);
  }
  
  void load_topodynamical_potential(topotential_pars_t& pars,
				    const bool& mandatory)
  {
    if(file_exists(topo_file_name)) pars.load(topo_file_name);
    else
      if(mandatory) crash("%s file not found when mandatory present",topo_file_name);
      else verbosity_lv2_master_printf("%s not found, skipping reading",topo_file_name);
  }
  
  //Compute the topological action
  double topotential_action(const EoField<quad_su3>& conf,
			    const topotential_pars_t &pars)
  {
    crash("reimplent");
    
    // //compute topocharge
    // double Q;
    // if(pars.stout_pars.nlevels)
    //   {
    // 	EoField<quad_su3> smeConf("smeConf",WITH_HALO);
    //     stout_smear(smeConf,conf,pars.stout_pars);
    // 	total_topological_charge_eo_conf(&Q,smeConf);
    //   }
    // else
    //   total_topological_charge_eo_conf(&Q,conf);
    
    //compute according to flag
    double topo_action=0;
    // switch(pars.flag)
    //   {
    //   case 1: topo_action=Q*pars.theta;break;
    //   case 2: topo_action=topodynamical_potential(Q,pars);break;
    //   default: crash("unknown flag %d",pars.flag);
    //   }
    
    // //free if it was allocated
    // if(pars.stout_pars.nlevels!=0) for(int eo=0;eo<2;eo++) nissa_free(conf[eo]);
    
    return topo_action;
  }
  
  //lx version
  double topotential_action(const LxField<quad_su3>& lx_conf,
			    const topotential_pars_t &pars)
  {
    //allocate
    EoField<quad_su3> eo_conf("stout_conf",WITH_HALO_EDGES);
    
    //split and compute
    split_lx_vector_into_eo_parts(eo_conf,lx_conf);
    
    return topotential_action(eo_conf,pars);
  }
  
  std::string topotential_pars_t::get_str(const bool& full) const
  {
    std::ostringstream os;
    
    const char name_known[3][10]={"NONE","ORDINARY","META"};
    
    if(full||flag!=def_flag()) os<<"TopoPotential\t=\t"<<name_known[flag]<<"\n";
    
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
}
