#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <base/export_conf_to_external_lib.hpp>
#include <base/multiGridParams.hpp>

namespace nissa
{
  bool export_gauge_conf_to_external_lib(quad_su3 *conf)
  {
    static checksum check_old={0,0},check_cur;
    
    //compute checksum
    checksum_compute_nissa_data(check_cur,conf,sizeof(quad_su3),sizeof(double)*8);
    
    //verify if export needed
    bool export_needed=false;
    for(int i=0;i<2;i++)
      {
	//check inited
	bool export_since_new=(check_old[i]==0);
	if(not export_needed and export_since_new) master_printf("external library: Old checksum 0, need to export the conf\n");
	export_needed|=export_since_new;
	//check diff
	bool export_since_diff=(check_old[i]!=check_cur[i]);
	if(not export_needed and export_since_diff) master_printf("external library Old checksum %d is %x, new is %x, need to import\n",i,check_old[i],check_cur[i]);
	export_needed|=export_since_diff;
	//save
	check_old[i]=check_cur[i];
      }
    
//     if(export_needed)
//       {
// 	bool export_result;
// 	double plaq=0.0;
	
// #ifdef USE_DDALPHAAMG
// 	DDalphaAMG_set_configuration((double*)conf,&DD::status);
// 	export_result=not DD::status.success;
// #endif
	
#ifdef USE_QUDA
// 	export_result=true;
	// plaq=
	  quda_iface::load_conf(conf);
#endif
	
// 	multiGrid::setup_valid=false;
	
// 	if(export_result)
// 	  verbosity_lv1_master_printf("external library conf set, plaquette %e\n",plaq);
// 	else
// 	  crash("configuration updating did not run correctly");
//       }
//     else master_printf("No import needed\n");
    
    return export_needed;
  }
}
