#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <algorithm>

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "io/input.hpp"
#include "new_types/su3.hpp"
#include "operations/smearing/stout.hpp"
#include "operations/su3_paths/topological_charge.hpp"
#include "routines/ios.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

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
}
