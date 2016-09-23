#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_BRIDGE
#include "tmLQCD_bridge.hpp"

#include "routines/ios.hpp"

namespace nissa
{
  //initialize tmLQCD
  void tmLQCD_init()
  {
    //print message
    verbosity_lv2_master_printf("Calling tmLQCD\n");
    
    //parameters
    int argc=0,verbose=1,external_id=0;
    char **argv=NULL;
    
    //call
    tmLQCD::tmLQCD_invert_init(argc,argv,verbose,external_id);
  }
  
  //export a gauge configuration to tmLQCD
  void export_gauge_conf_to_tmLQCD(quad_su3 *conf_lx)
  {
    external_conf_to_tmLQCD_handle=conf_lx;
    tmLQCD_import_gauge(nissa_feed_conf_to_tmLQCD);
  }
}
