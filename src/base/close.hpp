#ifndef _CLOSE_HPP
#define _CLOSE_HPP

#include <base/memory_manager.hpp>
#include <base/random.hpp>
#include <base/vectors.hpp>
#include <geometry/geometry_eo.hpp>
#include <geometry/geometry_lx.hpp>
#include <operations/remap_vector.hpp>

#ifdef USE_QUDA
# include <base/quda_bridge.hpp>
#endif

namespace nissa
{
  inline void close_nissa()
  {
    master_printf("Closing nissa\n");
    
    //unset remappers
    for(int mu=0;mu<NDIM;mu++)
      {
	if(remap_lx_to_locd[mu]) delete remap_lx_to_locd[mu];
	if(remap_locd_to_lx[mu]) delete remap_locd_to_lx[mu];
      }
    
#ifdef USE_QUDA
    if(use_quda) quda_iface::finalize();
#endif
     
    //unset lx geometry
    if(lxGeomInited) unset_lx_geometry();
    
    //unset eo geometry
    if(eo_geom_inited) unset_eo_geometry();
    
    //stop the random generator
    if(loc_rnd_gen_inited) stop_loc_rnd_gen();
    
    //print information over the maximum amount of memory used
    master_printf("Maximal memory used during the run: %zu bytes (",max_required_memory);
    if(is_master_rank()) fprintf_friendly_filesize(stdout,max_required_memory);
    master_printf(") per rank\n\n");
    
    //check wether there are still allocated vectors
    if(main_vect.next!=NULL && is_master_rank())
      {
	printf("Warning, there are still allocated vectors:\n");
	print_all_vect_content();
	printf("For a total of %zu bytes\n",compute_vect_memory_usage());
      }
    
    delete cpuMemoryManager;
#ifdef USE_CUDA
    delete gpuMemoryManager;
#endif
    
    tot_time+=take_time();
    master_printf("Total time: %lg s\n",tot_time);
#ifdef COMM_BENCH
    master_printf("Total communication time: %lg s\n",tot_comm_time);
#endif
    
#ifdef USE_CUDA
    storeTunedKernelInfo();
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    master_printf("   Ciao!\n\n");
    MPI_Finalize();
  }
}

#endif
