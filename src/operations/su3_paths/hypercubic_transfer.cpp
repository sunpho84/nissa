#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //return the path transferring the origin of the hypercube (c_hyp_ori) to the passed point 
  //(c_base+c_hyp_) (so that all link are daggered)
  void get_covariant_transport_to_hypercube_origin(su3 path,coords c_hyp_ori,const coords c_hyp_red,quad_su3 **conf)
  {
    //reset the patht to zero
    int nsubpaths=0;
    su3_put_to_zero(path);
    
    for(int mu=0;mu<4;mu++)
      if(c_hyp_red[mu]==1)
	{
	  //increment the number of subpath added
	  nsubpaths++;
	  
	  //move one step
	  coords cprime_hyp_red={c_hyp_red[0],c_hyp_red[1],c_hyp_red[2],c_hyp_red[3]};
	  cprime_hyp_red[mu]--;
	  
	  //take the full coordinates of moved point
	  coords cprime;
	  for(int nu=0;nu<4;nu++) cprime[nu]=c_hyp_ori[nu]+cprime_hyp_red[nu];
	  
	  //take lx index and convert it to eo
	  int cprime_lx=loclx_of_coord(cprime);
	  int cprime_eo=loceo_of_loclx[cprime_lx];
	  int cprime_parity=loclx_parity[cprime_lx];
	  
	  //take sub-path, attach proper link and summ it
	  su3 subpath;
	  get_covariant_transport_to_hypercube_origin(subpath,c_hyp_ori,cprime_hyp_red,conf);
	  su3_summ_the_dag_prod_su3(path,conf[cprime_parity][cprime_eo][mu],subpath);
	}
    
    //average or put to id if no path really involved
    if(nsubpaths) su3_prod_double(path,path,1.0/nsubpaths);
    else          su3_put_to_id(path);
  }
  void get_covariant_transport_to_hypercube_origin(su3 path,int ivol,int hyp_red,quad_su3 **conf)
  {
    coords c_hyp_red;
    red_coords_of_hypercubic_red_point(c_hyp_red,hyp_red);
    
    get_covariant_transport_to_hypercube_origin(path,loc_coord_of_loclx[ivol],c_hyp_red,conf);
  }
  
  //covariantly transfer the origin of each hypercube to the site marked as hyp_red
  void gauge_transfer_in_hypercube_from_origin(color **out,quad_su3 **conf,int hyp_red,color **in)
  {
    GET_THREAD_ID();
    
    //check that out is not in, and reset it
    if(out==in) crash("not possible");
    vector_reset(out[EVN]);
    vector_reset(out[ODD]);
    
    //get hypercubic reduced coords
    coords c_hyp_red;
    red_coords_of_hypercubic_red_point(c_hyp_red,hyp_red);
    
    //loop over all 2^4 hypercubes
    NISSA_PARALLEL_LOOP(hyp_cube,0,loc_vol/16)
      {
	//take vertex coords and lx
	coords c_hyp_ori;
	lx_coords_of_hypercube_vertex(c_hyp_ori,hyp_cube);
	
	//get corresponding coords of destination point
	coords c;
	for(int mu=0;mu<4;mu++) c[mu]=c_hyp_ori[mu]+c_hyp_red[mu];
	
	//find the destination and convert to eo
	int clx=loclx_of_coord(c);
	int ceo=loceo_of_loclx[clx];
	int cpar=loclx_parity[clx];
	
	//get the parallel transport and apply it
	su3 path;
	get_covariant_transport_to_hypercube_origin(path,c_hyp_ori,c_hyp_red,conf);
	unsafe_su3_prod_color(out[cpar][ceo],path,in[EVN][loceo_of_loclx[loclx_of_coord(c_hyp_ori)]]);
      }
    
    set_borders_invalid(out[EVN]);
    set_borders_invalid(out[ODD]);
  }
}
