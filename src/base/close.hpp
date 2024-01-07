#ifndef _CLOSE_HPP
#define _CLOSE_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <base/lattice.hpp>
#include <base/memory_manager.hpp>
#include <communicate/communicate.hpp>
#include <routines/mpiRoutines.hpp>

#ifdef USE_QUDA
# include <base/quda_bridge.hpp>
#endif

namespace nissa
{
  inline void closeNissa()
  {
    masterPrintf("Closing nissa\n");
    
#ifdef USE_QUDA
    if(use_quda) quda_iface::finalize();
#endif
    
    freeCommunicationBuffers();
    
    lat.reset();
    delete _lat;
    
    delete cpuMemoryManager;
#ifdef USE_CUDA
    delete gpuMemoryManager;
#endif
    
    mpiRanksBarrier();
    masterPrintf("   Ciao!\n\n");
    mpiFinalize();
  }
}

#endif
