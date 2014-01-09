#include "gauge_sweep.hpp"

#include "base/global_variables.hpp"
#include "base/debug.hpp"
#include "geometry/geometry_lx.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{  
  //constructor
  gauge_sweep_t::gauge_sweep_t()
  {
    //allocate space to convert from the geometry to be used (to be externalized)
    ivol_of_box_dir_par=nissa_malloc("ivol_of_box_dir_par",4*loc_vol,int);
    
    //mark not to have inited geom and paths
    comm_init_time=0;
    path_inited=par_geom_inited=false;
    max_cached_link=max_sending_link;
  }
  
  //destructor
  gauge_sweep_t::~gauge_sweep_t()
  {
    //for sure (apart from when it will be externalized) we must deallocate geometry constructor
    nissa_free(ilink_per_paths);
    
    //then if we inited geometry really, also site counters
    if(par_geom_inited) nissa_free(nsite_per_box_dir_par);
    
    //and path computer
    if(path_inited)
      {
	nissa_free(buf_out);
	nissa_free(buf_in);
	nissa_free(ivol_of_box_dir_par);
	for(int ibox=0;ibox<16;ibox++) delete box_comm[ibox];
      }
  }
  
  //initialize the geometry of the boxes subdirpar sets
  void gauge_sweep_t::init_box_dir_par_geometry(int ext_gpar,int(*par_comp)(coords ivol_coord,int dir))
  {
    comm_init_time=-take_time();
    
    //mark as inited, store parity and allocate each subbox numberer
    gpar=ext_gpar;
    par_geom_inited=true;
    nsite_per_box_dir_par=nissa_malloc("nsite_per_box_dir_par",16*4*gpar,int);
    
    //find the elements of all boxes
    int ibox_dir_par=0; //order in the parity-split order
    for(int ibox=0;ibox<16;ibox++)
      for(int dir=0;dir<4;dir++)
	for(int par=0;par<gpar;par++)
	  {
	    nsite_per_box_dir_par[par+gpar*(dir+4*ibox)]=0;
	    for(int isub=0;isub<nsite_per_box[ibox];isub++)
	      {
		//get coordinates of site
		coords isub_coord;
		coord_of_lx(isub_coord,isub,box_size[ibox]);
		
		//get coords in the local size, and parity
		coords ivol_coord;
		for(int mu=0;mu<4;mu++) ivol_coord[mu]=box_size[0][mu]*box_coord[ibox][mu]+isub_coord[mu];
		int site_par=par_comp(ivol_coord,dir);
		if(site_par>=gpar||site_par<0) crash("obtained par %d while expecting in the range [0,%d]",par,gpar-1);
		
		//map sites in current parity
		if(site_par==par)
		  {
		    ivol_of_box_dir_par[ibox_dir_par++]=lx_of_coord(ivol_coord,loc_size);
		    nsite_per_box_dir_par[par+gpar*(dir+4*ibox)]++;
		  }
	      }
	  }  
  }
  
  //init the list of paths
  void gauge_sweep_t::init_paths(int ext_nlinks_per_paths_site,void(*ext_add_paths_per_site_dir)(int *ilink_to_be_used,all_to_all_gathering_list_t &gat,int ivol,int mu))
  {
    //take external nlinks and mark
    nlinks_per_paths_site=ext_nlinks_per_paths_site;
    add_paths_per_site_dir=ext_add_paths_per_site_dir;
    path_inited=true;
    
    //allocate
    ilink_per_paths=nissa_malloc("ilink_per_paths",nlinks_per_paths_site*4*loc_vol,int);
    all_to_all_gathering_list_t *gl[16];
    for(int ibox=0;ibox<16;ibox++) gl[ibox]=new all_to_all_gathering_list_t;
    add_paths(gl);
    
    //initialize the communicator
    for(int ibox=0;ibox<16;ibox++)
      {
	box_comm[ibox]=new all_to_all_comm_t(*(gl[ibox]));
	delete gl[ibox];
      }
    
    //compute the maximum number of link to send and receive and allocate buffers
    for(int ibox=0;ibox<16;ibox++)
      {
	max_cached_link=std::max(max_cached_link,box_comm[ibox]->nel_in);
	max_sending_link=std::max(max_sending_link,box_comm[ibox]->nel_out);
      }
    buf_out=nissa_malloc("buf_out",max_sending_link,su3);
    buf_in=nissa_malloc("buf_in",max_cached_link,su3);
    
    //check cached
    verbosity_lv3_master_printf("Max cached links: %d\n",max_cached_link);
    if(max_cached_link>bord_vol+edge_vol) crash("larger buffer needed");
  }
  
  //add all the links needed to compute staple separately for each box
  THREADABLE_FUNCTION_2ARG(add_paths_to_gauge_sweep, gauge_sweep_t*,gs, all_to_all_gathering_list_t**,gl)
  {
    GET_THREAD_ID();
    
    //add the links to paths
    NISSA_PARALLEL_LOOP(ibox,0,16)
      {
	//find base for curr box
	int ibase=0;
	for(int jbox=0;jbox<ibox;jbox++) ibase+=nsite_per_box[jbox];
	ibase*=4;
	
	//scan all the elements of sub-box, selecting only those with the good parity
	for(int dir=0;dir<4;dir++)
	  for(int par=0;par<gs->gpar;par++)
	    {	  
	      for(int ibox_dir_par=ibase;ibox_dir_par<ibase+gs->nsite_per_box_dir_par[par+gs->gpar*(dir+4*ibox)];
		  ibox_dir_par++)
		{
		  int ivol=gs->ivol_of_box_dir_par[ibox_dir_par];
		  gs->add_paths_per_site_dir(gs->ilink_per_paths+ibox_dir_par*gs->nlinks_per_paths_site,*(gl[ibox]),ivol,dir);
		}
	      ibase+=gs->nsite_per_box_dir_par[par+gs->gpar*(dir+4*ibox)];
	    }
    }
    THREAD_BARRIER();
  }}

  //wrapper to use threads
  void gauge_sweep_t::add_paths(all_to_all_gathering_list_t **gl)
  {
    add_paths_to_gauge_sweep(this,gl);
  }
}
