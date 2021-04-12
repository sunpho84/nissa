#include "nissa.hpp"

namespace nissa
{
  //multiply the whole conf for stag phases
  void addrem_stagphases_to_lx_conf(quad_su3* lx_conf)
  {
    //we must ensure that nobody is using the conf
    THREAD_BARRIER();
    
    //work also on borders and edges if allocated and valid
    LocLxSite ending=locVol;
    if(check_borders_allocated(lx_conf,0) and check_borders_valid(lx_conf)) ending=locVolWithBord;
    if(check_edges_allocated(lx_conf,0) and check_edges_valid(lx_conf)) ending=locVolWithBordAndEdge;
    
    NISSA_PARALLEL_LOOP(ivol,0,ending)
      {
	int d=0;
	
	//phase in direction 1 is always 0 so nothing has to be done in that dir
	//if(d%2==1) su3_prod_double(lx_conf[ivol.nastyConvert()][1],lx_conf[ivol.nastyConvert()][1],-1);
	
	//direction 2
	d+=glbCoordOfLoclx[ivol.nastyConvert()][1];
	if(d%2==1) su3_prod_double(lx_conf[ivol.nastyConvert()][2],lx_conf[ivol.nastyConvert()][2],-1);
	
	//direction 3
	d+=glbCoordOfLoclx[ivol.nastyConvert()][2];
	if(d%2==1) su3_prod_double(lx_conf[ivol.nastyConvert()][3],lx_conf[ivol.nastyConvert()][3],-1);
	
	//direction 0
	d+=glbCoordOfLoclx[ivol.nastyConvert()][3];
	
	//putting the anti-periodic condition on the temporal border
	if(glbCoordOfLoclx[ivol.nastyConvert()][0]==glbSize[0]-1) d+=1;
	if(d%2==1) su3_prod_double(lx_conf[ivol.nastyConvert()][0],lx_conf[ivol.nastyConvert()][0],-1);
      }
    NISSA_PARALLEL_LOOP_END;
    
    THREAD_BARRIER();
  }
}
