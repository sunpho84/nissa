#ifndef _CLOSE_HPP
#define _CLOSE_HPP

#include <base/memory_manager.hpp>
#include <routines/mpiRoutines.hpp>

#ifdef USE_QUDA
# include <base/quda_bridge.hpp>
#endif

namespace nissa
{
  inline void closeNissa()
  {
    master_printf("Closing nissa\n");
    
#ifdef USE_QUDA
    if(use_quda) quda_iface::finalize();
#endif
    
    //unset lx geometry
    CRASH("lat to be deleted");
    
    delete cpuMemoryManager;
#ifdef USE_CUDA
    delete gpuMemoryManager;
#endif
    
    mpiRanksBarrier();
    master_printf("   Ciao!\n\n");
    mpiFinalize();
  }
}

#endif
