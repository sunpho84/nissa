#ifndef _ARBITRARY_HPP
#define _ARBITRARY_HPP

#include <stdlib.h>
#include <deque>

#include <base/debug.hpp>
#include <base/vectors.hpp>
#include <geometry/geometry_lx.hpp>
#include <routines/ios.hpp>

#define START_PATH_FLAG 1
#define DAG_LINK_FLAG 2
#define NONLOC_LINK_FLAG 4
#define SUMM_TO_PREVIOUS_PATH_FLAG 8
#define nposs_path_flags 4

namespace nissa
{
  struct movement_link_id
  {
    int mov;
    int link_id;
    int ord;
    void set(int ext_mov,int gx,int mu){
      mov=ext_mov;
      link_id=gx*4+mu;
      int lx,rx;
      get_loclx_and_rank_of_glblx(&lx,&rx,gx);
      ord=((rank+nranks-rx)%nranks*locVol.nastyConvert()+lx)*4+mu; //sort according to recv rank
    }
  };
  
  struct paths_calculation_structure
  {
    //initialization
    paths_calculation_structure(int ext_npaths,int ext_ntot_mov) {
      //set npaths and tot mov, start position and move
      npaths=ext_npaths;
      ntot_mov=ext_ntot_mov;
      cur_path=0;
      cur_mov=0;
      
      verbosity_lv3_master_printf("Initializing a new path calculation structure with %d movements and %d paths\n",ntot_mov,npaths);
      link_for_movements=nissa_malloc("link_for_movements",ntot_mov,int);
      
      //mark that we have not finished last path
      finished_last_path_flag=0;
      
      //we do not yet know hown many ranks we need to communicate with
      nranks_to_send=nranks_to_recv=0;
      nnonloc_links=nind_nonloc_links=0;
      
      //and we don't know which they are
      ranks_to_send_list=ranks_to_recv_list=NULL;
      
      //reset the list of links to send, and the amount of links to send and receive from each rank
      links_to_send_list=NULL;
      nlinks_to_recv_list=nlinks_to_send_list=NULL;
      
      //reset the list of pairs of link and movements
      movements_nonloc_links_id_list=NULL;
    }
    ~paths_calculation_structure() {
      nissa_free(link_for_movements);
      
      if(finished_last_path_flag)
	{
	  nissa_free(ind_nonloc_links_list);
	  
	  //free the list of links to receive
	  free(links_to_send_list);
	  
	  //the amount of links to send and receive from each rank
	  free(nlinks_to_send_list);
	  free(nlinks_to_recv_list);
	  
	  //and the list of ranks to communicate with
	  free(ranks_to_send_list);
	  free(ranks_to_recv_list);
	}
    }
    
    //parameters defining the set of paths
    int npaths,ntot_mov;
    int *link_for_movements;
    movement_link_id *movements_nonloc_links_id_list;
    
    //current global movement (link), path and position
    int cur_path,cur_mov;
    int pos;
    
    //relevant for MPI part
    int finished_last_path_flag;
    int nnonloc_links,nind_nonloc_links;
    int nranks_to_send,*ranks_to_send_list,*nlinks_to_send_list,ntot_links_to_send;
    int nranks_to_recv,*ranks_to_recv_list,*nlinks_to_recv_list,ntot_links_to_recv;
    int *links_to_send_list,*ind_nonloc_links_list;
    
    //commands
    void start_new_path_from_loclx(const LocLxSite& lx)
    {
      pos=glblxOfLoclx(lx).nastyConvert();
      link_for_movements[cur_mov]=START_PATH_FLAG;
    };
    
    void summ_to_previous_path()
    {
      if(cur_path==0)
	crash("cannot summ to path number 0");
      
      //if not already set, diminuish the number of paths
      if(!(link_for_movements[cur_mov]&SUMM_TO_PREVIOUS_PATH_FLAG)) cur_path--;
      link_for_movements[cur_mov]+=SUMM_TO_PREVIOUS_PATH_FLAG;
    }
    
    void switch_to_next_step()
    {
      cur_mov++;
      if(cur_mov>ntot_mov) crash("exceded (%d) the number of allocatec movements, %d",cur_mov,ntot_mov);
      link_for_movements[cur_mov]=0;
    }
    
    void move_forward(int mu);
    void move_backward(int mu);
    void stop_current_path()
    {
      cur_path++;
      if(cur_path>npaths) crash("exceded (%d) the number of allocated paths, %d",cur_path,npaths);
    }
    
    void finished_last_path();
    void gather_nonloc_start(MPI_Request *request,int &irequest,su3 *nonloc_links);
    void gather_nonloc_finish(MPI_Request *request,int &irequest,su3 *send_buff);
    void gather_nonloc_lx(su3 *paths,quad_su3 *conf);
    void gather_nonloc_eo(su3 *paths,quad_su3 **conf);
    void compute_lx(su3 *paths,quad_su3 *conf);
    void compute_eo(su3 *paths,quad_su3 **conf);
    
  private:
    //avoid bare initialization without specification of nel
    paths_calculation_structure();
    //setu the communication buffers
    void setup_sender_receivers();
  };
  
  ////////////////////////////////////////////////////////////
  
  struct coords_t{
    coords c;
    int &operator[](int i){return c[i];}
    bool operator==(coords_t in){bool out=true;for(int i=0;i<NDIM;i++) out&=(c[i]==in.c[i]);return out;}
    bool operator!=(coords_t i){return !((*this)==i);}
    coords_t(){memset(c,0,sizeof(coords));}
    coords_t(const coords_t &o){memcpy(c,o.c,sizeof(coords));}
  };
  typedef std::deque<coords_t> path_drawing_t;
  
  void init_su3_path(path_drawing_t *c,su3 *out);
  void elong_su3_path_BW(path_drawing_t *c,su3 *out,quad_su3 *conf,int mu,bool both_sides=false);
  void elong_su3_path_FW(path_drawing_t *c,su3 *out,quad_su3 *conf,int mu,bool both_sides=false);
  void elong_su3_path(path_drawing_t *c,su3 *out,quad_su3 *conf,int mu,int len,bool both_sides=false);
  
  //direction and length
  typedef std::pair<int,int> path_step_pars_t;
  //list of direction,lengths
  typedef std::vector<path_step_pars_t> path_list_steps_t;
  
  //wrapper
  inline void elong_su3_path(path_drawing_t *c,su3 *out,quad_su3 *conf,path_step_pars_t pars)
  {elong_su3_path(c,out,conf,pars.first,pars.second);}
  void elong_su3_path(path_drawing_t *c,su3 *out,quad_su3 *conf,path_list_steps_t steps);
  inline void compute_su3_path(path_drawing_t *c,su3 *out,quad_su3 *conf,path_list_steps_t steps)
  {
    init_su3_path(c,out);
    elong_su3_path(c,out,conf,steps);
    
  }
}

#endif

