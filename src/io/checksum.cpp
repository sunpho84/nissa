#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <stdint.h>
#include <unistd.h>

#include "checksum.hpp"

#include "base/debug.hpp"
#include "geometry/geometry_lx.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //compute the checksum of ildg data (little endianness, time is slower index)
  void checksum_compute_ildg_data(uint32_t *check,void *data,size_t bps)
  {
    const double init_time=take_time();
    
    uint32_t loc_check[2]={0,0};
    
    NISSA_LOC_VOL_LOOP(ivol)
      {
	const coords_t& x=locCoordOfLoclx[ivol];
	const coords_t& X=glbCoordOfLoclx[ivol];
	uint32_t loc_ivol=x[0],glb_ivol=X[0];
	for(int mu=NDIM-1;mu>0;mu--)
	  {
	    loc_ivol=loc_ivol*locSize[mu]+x[mu];
	    glb_ivol=glb_ivol*glbSize[mu]+X[mu];
	  }
	uint32_t crc_rank[2]={glb_ivol%29,glb_ivol%31};
	
	uint32_t temp=ildg_crc32(0,(unsigned char*)data+bps*loc_ivol,bps);
	
	for(int i=0;i<2;i++) loc_check[i]^=temp<<crc_rank[i]|temp>>(32-crc_rank[i]);
      }
    
    MPI_Allreduce(loc_check,check,2,MPI_UNSIGNED,MPI_BXOR,MPI_COMM_WORLD);
    
    master_printf("time to compute checksum using cpp version: %lg (%s), %zu bps\n",take_time()-init_time,__PRETTY_FUNCTION__,bps);
  }
}
