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
    enum ExportBypass{NO_BYPASS,AVOID_EXPORT,FORCE_EXPORT};
    EXTERN_EXPORT_CONF ExportBypass export_bypass INIT_EXPORT_CONF_TO(=NO_BYPASS);
    EXTERN_EXPORT_CONF checksum check_old INIT_EXPORT_CONF_TO(={0,0});
    EXTERN_EXPORT_CONF std::string confTagOld INIT_EXPORT_CONF_TO(="");
    EXTERN_EXPORT_CONF std::string confTag INIT_EXPORT_CONF_TO(="");
  }
  
  /// Keep track of the exported conf
  template <typename T>
  bool export_gauge_conf_to_external_solver(const T& conf)
  {
    using namespace export_conf;
    
    if(export_bypass==AVOID_EXPORT)
      return false;
    
    //verify if export needed
    bool export_needed=false;
    
    if(confTag!=confTagOld)
      {
	master_printf("previously exported conf tag: \"%s\", to be exported: \"%s\", bypassing the checksum\n",confTagOld.c_str(),confTag.c_str());
	
	export_conf::export_bypass=
	  export_conf::FORCE_EXPORT;
	
	confTagOld=confTag;
      }
    
    if(export_bypass==FORCE_EXPORT)
      {
	master_printf("Forcing export of the conf to external library\n");
	export_needed=true;
	export_bypass=NO_BYPASS;
      }
    else
	{
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
