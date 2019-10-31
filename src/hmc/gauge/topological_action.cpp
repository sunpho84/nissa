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
  double topodynamical_potential(double Q,topotential_pars_t &pars)
  {return pars.compute_pot(Q);}
  //draw the topodynamical potential
  void save_topodynamical_potential(topotential_pars_t &pars)
  {pars.save(topo_file_name);}
  void load_topodynamical_potential(topotential_pars_t &pars,bool mandatory)
  {
    if(file_exists(topo_file_name)) pars.load(topo_file_name);
    else
      if(mandatory) crash("%s file not found when mandatory present",topo_file_name);
      else verbosity_lv2_master_printf("%s not found, skipping reading",topo_file_name);
  }
  
  //Compute the topological action
  double topotential_action(quad_su3 **ext_conf,topotential_pars_t &pars)
  {
    quad_su3 *conf[2];
    if(pars.stout_pars.nlevels==0)
      {
        conf[0]=ext_conf[0];
        conf[1]=ext_conf[1];
      }
    else
      {
        conf[0]=nissa_malloc("stout_conf",loc_volh+bord_volh+edge_volh,quad_su3);
        conf[1]=nissa_malloc("stout_conf",loc_volh+bord_volh+edge_volh,quad_su3);
        
        //smear
        stout_smear(conf,ext_conf,&(pars.stout_pars));
      }
    
    //compute topocharge
    double Q;
    total_topological_charge_eo_conf(&Q,conf);
    
    //compute according to flag
    double topo_action=0;
    switch(pars.flag)
      {
      case 1: topo_action=Q*pars.theta;break;
      case 2: topo_action=topodynamical_potential(Q,pars);break;
      default: crash("unknown flag %d",pars.flag);
      }
    
    //free if it was allocated
    if(pars.stout_pars.nlevels!=0) for(int eo=0;eo<2;eo++) nissa_free(conf[eo]);
    
    return topo_action;
  }
  
  //lx version
  double topotential_action(quad_su3 *lx_conf,topotential_pars_t &pars)
  {
    //allocate
    quad_su3 *eo_conf[2];
    for(int eo=0;eo<2;eo++) eo_conf[eo]=nissa_malloc("stout_conf",loc_volh+bord_volh+edge_volh,quad_su3);
    
    //split and compute
    split_lx_vector_into_eo_parts(eo_conf,lx_conf);
    double out=topotential_action(eo_conf,pars);
    
    //free and return
    for(int eo=0;eo<2;eo++) nissa_free(eo_conf[eo]);
    return out;
  }
  
  std::string topotential_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    const char name_known[3][10]={"NONE","ORDINARY","META"};
    if(full||flag!=def_flag()) os<<"TopoPotential\t=\t"<<name_known[flag]<<"\n";
    switch(flag)
      {
      case 0:break;
      case 1:os<<" Theta\t\t=\t"<<theta<<"\n";break;
      case 2:
	os<<meta_pars_t::get_str(full);
	os<<stout_pars.get_str(full);
	break;
      }
    return os.str();
    }
}
