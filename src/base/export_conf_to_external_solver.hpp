#ifndef _EXPORT_CONF_TO_EXTERNAL_SOLVER_HPP
#define _EXPORT_CONF_TO_EXTERNAL_SOLVER_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef USE_DDALPHAAMG
# include "base/DDalphaAMG_bridge.hpp"
#endif

#ifdef USE_QUDA
# include "base/quda_bridge.hpp"
#endif

#include "io/checksum.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"

#ifndef EXTERN_EXPORT_CONF
 #define EXTERN_EXPORT_CONF extern
 #define INIT_EXPORT_CONF_TO(cond...)
#else
 #define INIT_EXPORT_CONF_TO(cond...) cond
#endif

namespace nissa
{
  namespace export_conf
  {
    enum ExportRule{DO_THE_CHECK,FORCE_EXPORT,AVOID_EXPORT};
    EXTERN_EXPORT_CONF ExportRule export_rule INIT_EXPORT_CONF_TO(=DO_THE_CHECK);
    EXTERN_EXPORT_CONF bool relyOnTag INIT_EXPORT_CONF_TO(=false);
    EXTERN_EXPORT_CONF checksum check_old INIT_EXPORT_CONF_TO(={0,0});
    EXTERN_EXPORT_CONF std::string confTagOld INIT_EXPORT_CONF_TO(="");
    EXTERN_EXPORT_CONF std::string confTag INIT_EXPORT_CONF_TO(="");
  }
  
  /// Keep track of the exported conf
  template <typename T>
  bool export_gauge_conf_to_external_solver(const T& conf)
  {
    using namespace export_conf;
    
    if(export_rule==AVOID_EXPORT)
      return false;
    
    //verify if export needed
    bool export_needed=false;
    if(export_rule==FORCE_EXPORT)
      {
	master_printf("Forcing export of the conf to external library\n");
	export_needed=true;
      }
    else
      if(relyOnTag)
	{
	  master_printf("Relying on tag to check,\n old tag: \"%s\"\n",confTagOld.c_str());
	  
	  if(confTag!=confTagOld)
	    {
	      master_printf("new tag: \"%s\"\n -> export needed\n",confTag.c_str());
	      
	      confTagOld=confTag;
	      export_needed=true;
	    }
	  else
	    master_printf(" -> tag not changed, avoiding export\n",confTag.c_str());
	}
      else
	{
	  master_printf("Relying on checksum to check\n");
	  
	  checksum check_cur{};
	  
	  //compute checksum
	  checksum_compute_nissa_data(check_cur,conf,sizeof(double)*8,sizeof(quad_su3));
	  
	  for(int i=0;i<2;i++)
	    {
	      //check inited
	      bool export_since_new=(check_old[i]==0);
	      if(export_since_new) master_printf("external library: Old checksum 0, need to export the conf\n");
	      export_needed|=export_since_new;
	      //check diff
	      bool export_since_diff=(check_old[i]!=check_cur[i]);
	      if(export_since_diff) master_printf("external library Old checksum %d is %x, new is %x, need to import\n",i,check_old[i],check_cur[i]);
	      export_needed|=export_since_diff;
	      //save
	      check_old[i]=check_cur[i];
	    }
	}
    
    if(export_needed)
      {
	bool export_result=false;
	double plaq=0.0;
	
#ifdef USE_DDALPHAAMG
	DD::set_configuration(conf);
	export_result=DD::status.success;
#endif
	
#ifdef USE_QUDA
	export_result=true;
	plaq=
	  quda_iface::load_conf(conf);
#endif
	
	multiGrid::setup_valid=false;
	
	if(export_result)
	  verbosity_lv1_master_printf("external library conf set, plaquette %lg\n",plaq);
	else
	  crash("configuration updating did not run correctly");
      }
    else master_printf("No import needed\n");
    
    return export_needed;
  }
}

#undef INIT_EXPORT_CONF_TO
#undef EXTERN_EXPORT_CONF

#endif
