#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_TMLQCD_BRIDGE
#include "tmLQCD_bridge.hpp"

#include "routines/ios.hpp"
#include "routines/thread.hpp"

namespace nissa
{
  //initialize tmLQCD
  void tmLQCD_init()
  {
    //print message
    VERBOSITY_LV2_MASTER_PRINTF("Calling tmLQCD\n");
    
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
  
  //write the input file
  FILE* open_prepare_input_file_for_tmLQCD()
  {
    FILE *ftemp=open_file("invert.input","w");
    master_fprintf(ftemp,"L=%d\n",glb_size[1]);
    master_fprintf(ftemp,"T=%d\n",glb_size[0]);
    master_fprintf(ftemp,"NrXProcs=%d\n",nrank_dir[1]);
    master_fprintf(ftemp,"NrYProcs=%d\n",nrank_dir[2]);
    master_fprintf(ftemp,"NrZProcs=%d\n",nrank_dir[3]);
    master_fprintf(ftemp,"OMPNumThreads=%d\n",nthreads);
    
    // DisableIOChecks = yes
    // DebugLevel = 1
    // InitialStoreCounter = 1000
    // Measurements = 1
    // 2kappamu = 0.177
    // kappa = 0.177
    // GaugeConfigInputFile = conf
    // UseEvenOdd = yes
    
    return ftemp;
  }
}
