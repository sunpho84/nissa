#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <algorithm>

#include "base/thread_macros.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/new_types_definitions.hpp"
#include "new_types/su3.hpp"
#include "operations/smearing/stout.hpp"
#include "operations/su3_paths/topological_charge.hpp"
#include "routines/ios.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif


namespace nissa
{
  //compute the topodynamical potential using past history
  double topodynamical_potential(double Q,topotential_pars_t &pars)
  {
    double topotential=0;
    for(std::vector<double>::iterator it=pars.past_values.begin();it!=pars.past_values.end();it++)
      {
        double q=*it;
        double diff=Q-q,f=diff/pars.width,cont=exp(-f*f/2);
        topotential+=cont;
      }
    topotential*=pars.coeff;
    
    return topotential;
  }
  
  //draw the topodynamical potential
  void draw_topodynamical_potential(topotential_pars_t &pars)
  {
    FILE *file=open_file("topot","w");
    
    double Q_min=*std::min_element(pars.past_values.begin(),pars.past_values.end());
    double Q_max=*std::max_element(pars.past_values.begin(),pars.past_values.end());
    double Q_diff=Q_max-Q_min;
    int n=ceil(Q_diff/pars.width*20);
    double dQ=Q_diff/n;
    
    for(double Q=Q_min;Q<=Q_max;Q+=dQ) master_fprintf(file,"%lg %lg\n",Q,topodynamical_potential(Q,pars));
    close_file(file);
  }

  //Compute the topological action
  double topotential_action(quad_su3 **ext_conf,topotential_pars_t &pars)
  {
    quad_su3 *conf[2];
    if(pars.stout_pars.nlev==0)
      {
        conf[0]=ext_conf[0];
        conf[1]=ext_conf[1];
      }
    else
      {
        conf[0]=nissa_malloc("stout_conf",loc_volh+bord_volh+edge_volh,quad_su3);
        conf[1]=nissa_malloc("stout_conf",loc_volh+bord_volh+edge_volh,quad_su3);
        
        //smear
        addrem_stagphases_to_eo_conf(ext_conf);
        stout_smear(conf,ext_conf,&(pars.stout_pars));
        addrem_stagphases_to_eo_conf(ext_conf);
        addrem_stagphases_to_eo_conf(conf);
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
    if(pars.stout_pars.nlev!=0) for(int eo=0;eo<2;eo++) nissa_free(conf[eo]);
    
    return topo_action;
  }
}
