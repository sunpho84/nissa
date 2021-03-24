#include "nissa.hpp"

namespace nissa
{
  //multiply the whole conf for stag phases
  void addrem_stagphases_to_lx_conf(quad_su3* lx_conf)
  {
    //we must ensure that nobody is using the conf
    THREAD_BARRIER();
    
    //work also on borders and edges if allocated and valid
    int ending=loc_vol;
    if(check_borders_allocated(lx_conf,0) && check_borders_valid(lx_conf)) ending+=bord_vol;
    if(check_edges_allocated(lx_conf,0) && check_edges_valid(lx_conf)) ending+=edge_vol;
    
    NISSA_PARALLEL_LOOP(ivol,0,ending)
      {
	int d=0;
	
	//phase in direction 1 is always 0 so nothing has to be done in that dir
	//if(d%2==1) su3_prod_double(lx_conf[ivol][1],lx_conf[ivol][1],-1);
	
	//direction 2
	d+=glb_coord_of_loclx[ivol][1];
	if(d%2==1) su3_prod_double(lx_conf[ivol][2],lx_conf[ivol][2],-1);
	
	//direction 3
	d+=glb_coord_of_loclx[ivol][2];
	if(d%2==1) su3_prod_double(lx_conf[ivol][3],lx_conf[ivol][3],-1);
	
	//direction 0
	d+=glb_coord_of_loclx[ivol][3];
	
	//putting the anti-periodic condition on the temporal border
	if(glb_coord_of_loclx[ivol][0]==glb_size[0]-1) d+=1;
	if(d%2==1) su3_prod_double(lx_conf[ivol][0],lx_conf[ivol][0],-1);
      }
    NISSA_PARALLEL_LOOP_END;
    
    THREAD_BARRIER();
  }
}
