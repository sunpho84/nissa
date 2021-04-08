#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_mix.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3_op.hpp"
#include "operations/shift.hpp"
#include "routines/mpi_routines.hpp"
#include "threads/threads.hpp"

#include "arbitrary.hpp"

namespace nissa
{
  int compare_movement_link_id(const void *a,const void *b)
  {
    return ((movement_link_id*)a)->ord-((movement_link_id*)b)->ord;
  }
  
  //add a forward move
  void paths_calculation_structure::move_forward(int mu)
  {
    //check not to have passed the max number of steps
    if(cur_mov==ntot_mov) crash("exceded (%d) the number of allocated movements, %d",cur_mov,ntot_mov);
    //find global pos
    int nx=glblx_neighup(pos,mu);
    //search rank hosting site and loclx
    int lx,rx;
    get_loclx_and_rank_of_glblx(&lx,&rx,pos);
    //if local link, mark it, otherwise add to the list of non-locals
    if(rx==rank) link_for_movements[cur_mov]+=(lx*4+mu)<<nposs_path_flags;
    else
      {
	movements_nonloc_links_id_list=(movement_link_id*)realloc(movements_nonloc_links_id_list,sizeof(movement_link_id)*(nnonloc_links+1));
	movements_nonloc_links_id_list[nnonloc_links++].set(cur_mov,pos,mu);
      }
    //set new pos
    pos=nx;
    //switch to next step
    switch_to_next_step();
  }
  
  //add a backward move
  void paths_calculation_structure::move_backward(int mu)
  {
    //mark backward move
    link_for_movements[cur_mov]+=DAG_LINK_FLAG;
    //check not to have passed the max number of steps
    if(cur_mov==ntot_mov) crash("exceeded (%d) the number of allocated movements, %d",cur_mov,ntot_mov);
    //find global pos
    int nx=glblx_neighdw(pos,mu);
    //search rank hosting site and loclx
    int lx,rx;
    get_loclx_and_rank_of_glblx(&lx,&rx,nx);
    //if local link, mark it, otherwise add to the list of non-locals
    if(rx==rank) link_for_movements[cur_mov]+=(lx*4+mu)<<nposs_path_flags;
    else
      {
	movements_nonloc_links_id_list=(movement_link_id*)realloc(movements_nonloc_links_id_list,sizeof(movement_link_id)*(nnonloc_links+1));
	movements_nonloc_links_id_list[nnonloc_links++].set(cur_mov,nx,mu);
      }
    //set new pos
    pos=nx;
    //switch to next step
    switch_to_next_step();
  }
  
  //finish settings the paths, setup the send and receiver
  void paths_calculation_structure::finished_last_path()
  {
    if(cur_path!=npaths) crash("finished the path list at path %d while it was initialized for %d",cur_path,npaths);
    if(cur_mov!=ntot_mov) crash("finished the path list at mov %d while it was initialized for %d",cur_mov,ntot_mov);
    
    //sort the list
    qsort(movements_nonloc_links_id_list,nnonloc_links,sizeof(movement_link_id),compare_movement_link_id);
    
    //asign the non local links for each movement one by one, counting the independent ones
    nind_nonloc_links=0;
    for(int ilink=0;ilink<nnonloc_links;ilink++)
      {
	link_for_movements[movements_nonloc_links_id_list[ilink].mov]+=(nind_nonloc_links<<nposs_path_flags)+NONLOC_LINK_FLAG;
	if(ilink==(nnonloc_links-1)||movements_nonloc_links_id_list[ilink].link_id!=movements_nonloc_links_id_list[ilink+1].link_id) nind_nonloc_links++;
      }
    
    //allocate the list of nonlocal indep links to ask, so we can free the full list
    ind_nonloc_links_list=nissa_malloc("ind_nlonloc_links_list",nind_nonloc_links,int);
    nind_nonloc_links=0;
    for(int ilink=0;ilink<nnonloc_links;ilink++)
      {
	ind_nonloc_links_list[nind_nonloc_links]=movements_nonloc_links_id_list[ilink].link_id;
	if(ilink==(nnonloc_links-1)||movements_nonloc_links_id_list[ilink].link_id!=movements_nonloc_links_id_list[ilink+1].link_id)
	  nind_nonloc_links++;
      }
    
    //free here
    free(movements_nonloc_links_id_list);
    
    finished_last_path_flag=1;
    
    setup_sender_receivers();
  }
  
  //setup the sender and receiver buffers, finding which ranks are involved
  void paths_calculation_structure::setup_sender_receivers()
  {
    ntot_links_to_send=ntot_links_to_recv=0;
    
    //loop over the ranks displacement
    for(int delta_rank=1;delta_rank<nranks;delta_rank++)
      {
	int rank_to_send=(rank+delta_rank)%nranks;
	int rank_to_recv=(rank+nranks-delta_rank)%nranks;
	
	//counts the number of links to receive
	int nlinks_to_recv=0;
	for(int ilink=0;ilink<nind_nonloc_links;ilink++)
	  {
	    int t=ind_nonloc_links_list[ilink];
	    int gx=t>>2;
	    int rx=rank_hosting_glblx(gx);
	    if(rx==rank_to_recv) nlinks_to_recv++;
	  }
	
	//send this piece of info and receive the number of links to send
	int nlinks_to_send=0;
	MPI_Sendrecv((void*)&(nlinks_to_recv),1,MPI_INT,rank_to_recv,rank_to_recv*nranks+rank,
		     (void*)&(nlinks_to_send),1,MPI_INT,rank_to_send,rank*nranks+rank_to_send,
		     MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	
	//allocate a buffer where to store the list of links to ask
	int *links_to_ask=nissa_malloc("links_to_ask",nlinks_to_recv,int);
	nlinks_to_recv=0;
	for(int ilink=0;ilink<nind_nonloc_links;ilink++)
	  {
	    int t=ind_nonloc_links_list[ilink];
	    int gx=t>>2;
	    int mu=t%4;
	    
	    //get lx and rank hosting the site
	    int lx,rx;
	    get_loclx_and_rank_of_glblx(&lx,&rx,gx);
	    
	    //copy in the list if appropriate rank
	    if(rx==rank_to_recv)
	      {
		links_to_ask[nlinks_to_recv]=lx*4+mu;
		nlinks_to_recv++;
	      }
	  }
	
	//allocate the list of link to send
	int *links_to_send=nissa_malloc("links_to_send",nlinks_to_send,int);
	
	//send this piece of info
	MPI_Sendrecv((void*)links_to_ask, nlinks_to_recv,MPI_INT,rank_to_recv,rank_to_recv*nranks+rank,
		     (void*)links_to_send,nlinks_to_send,MPI_INT,rank_to_send,rank*nranks+rank_to_send,
		     MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	nissa_free(links_to_ask);
	
	//store the sending rank id and list of links
	if(nlinks_to_send!=0)
	  {
	    //prepare the space
	    ranks_to_send_list=(int*)realloc(ranks_to_send_list,sizeof(int)*(nranks_to_send+1));
	    links_to_send_list=(int*)realloc(links_to_send_list,sizeof(int)*(ntot_links_to_send+nlinks_to_send));
	    nlinks_to_send_list=(int*)realloc(nlinks_to_send_list,sizeof(int)*(nranks_to_send+1));
	    
	    //store
	    ranks_to_send_list[nranks_to_send]=rank_to_send;
	    memcpy(links_to_send_list+ntot_links_to_send,links_to_send,nlinks_to_send*sizeof(int));
	    nlinks_to_send_list[nranks_to_send]=nlinks_to_send;
	    
	    //increase
	    nranks_to_send++;
	    ntot_links_to_send+=nlinks_to_send;
	  }
	nissa_free(links_to_send);
	
	//store the receiving rank id and list of links
	if(nlinks_to_recv!=0)
	  {
	    //prepare the space
	    ranks_to_recv_list=(int*)realloc(ranks_to_recv_list,sizeof(int)*(nranks_to_recv+1));
	    nlinks_to_recv_list=(int*)realloc(nlinks_to_recv_list,sizeof(int)*(nranks_to_recv+1));
	    
	    //store
	    ranks_to_recv_list[nranks_to_recv]=rank_to_recv;
	    nlinks_to_recv_list[nranks_to_recv]=nlinks_to_recv;
	    
	    //increase
	    nranks_to_recv++;
	    ntot_links_to_recv+=nlinks_to_recv;
	  }
      }
  }
  
  //open incoming communications
  void paths_calculation_structure::gather_nonloc_start(MPI_Request *request,int &irequest,su3 *nonloc_links)
  {
    //open receiving communications
    su3 *recv_ptr=nonloc_links;
    for(int irecv=0;irecv<nranks_to_recv;irecv++)
      {
	MPI_Irecv((void*)recv_ptr,nlinks_to_recv_list[irecv],MPI_SU3,ranks_to_recv_list[irecv],rank*nranks+ranks_to_recv_list[irecv],cart_comm,request+irequest++);
	
	if(irecv+1!=nranks_to_recv) recv_ptr+=nlinks_to_recv_list[irecv];
      }
  }
  
  //communicate
  void paths_calculation_structure::gather_nonloc_finish(MPI_Request *request,int &irequest,su3 *send_buff)
  {
    //open sending communications
    su3 *send_ptr=send_buff;
    for(int isend=0;isend<nranks_to_send;isend++)
      {
	MPI_Isend((void*)send_ptr,nlinks_to_send_list[isend],MPI_SU3,ranks_to_send_list[isend],ranks_to_send_list[isend]*nranks+rank,cart_comm,request+irequest++);
	
	if(isend+1!=nranks_to_send) send_ptr+=nlinks_to_send_list[isend];
      }
    
    //wait communications to finish
    int ntot_req=nranks_to_send+nranks_to_recv;
    MPI_Waitall(ntot_req,request,MPI_STATUS_IGNORE);
  }
  
  //collect all the independent links entering into the calculations into "links"
  void paths_calculation_structure::gather_nonloc_lx(su3 *nonloc_links,quad_su3 *conf)
  {
    //list of request
    MPI_Request request[nranks_to_send+nranks_to_recv];
    int irequest=0;
    
    //open incoming communications
    gather_nonloc_start(request,irequest,nonloc_links);
    
    //allocate the sending buffer, fill them
    su3 *send_buff=nissa_malloc("buff",ntot_links_to_send,su3);
    for(int ilink=0;ilink<ntot_links_to_send;ilink++) su3_copy(send_buff[ilink],((su3*)conf)[links_to_send_list[ilink]]);
    
    //finish the communications
    gather_nonloc_finish(request,irequest,send_buff);
    
    nissa_free(send_buff);
  }
  
  //collect all the independent links entering into the calculations into "links" using the infamous eo geometry
  void paths_calculation_structure::gather_nonloc_eo(su3 *nonloc_links,quad_su3 **conf)
  {
    //list of request
    MPI_Request request[nranks_to_send+nranks_to_recv];
    int irequest=0;
    
    //open incoming communications
    gather_nonloc_start(request,irequest,nonloc_links);
    
    //allocate the sending buffer, fill them
    su3 *send_buff=nissa_malloc("buff",ntot_links_to_send,su3);
    for(int ilink=0;ilink<ntot_links_to_send;ilink++)
      {
	int t=links_to_send_list[ilink];
	int lx=t/4,mu=t%4;
	int p=loclx_parity[lx],eo=loceo_of_loclx[lx];
	su3_copy(send_buff[ilink],conf[p][eo][mu]);
      }
    
    //finish the communications
    gather_nonloc_finish(request,irequest,send_buff);
    
    nissa_free(send_buff);
  }
  
  void paths_calculation_structure::compute_lx(su3 *paths,quad_su3 *conf)
  {
    //buffer for non local links gathering
    su3 *nonloc_links=nissa_malloc("links",nind_nonloc_links,su3);
    
    //gather non local links
    gather_nonloc_lx(nonloc_links,conf);
    
    //running path
    su3 cur_path;
    su3_put_to_zero(cur_path);
    
    //compute the paths one by one
    int ipath=0;
    for(int imov=0;imov<ntot_mov;imov++)
      {
	int ilink=link_for_movements[imov]>>nposs_path_flags;
	
	//take the flags
	int tag=link_for_movements[imov]%(1<<nposs_path_flags);
	int start=tag&START_PATH_FLAG;
	int herm=tag&DAG_LINK_FLAG;
	int nonloc=tag&NONLOC_LINK_FLAG;
	int summ_to_previous=tag&SUMM_TO_PREVIOUS_PATH_FLAG;
	
	//check if starting a new path
	if(start==1)
	  {
	    //if not the first mov, start the new path
	    if(imov!=0)
	      {
		//summ the ultimated path in its place
		su3_summassign(paths[ipath],cur_path);
		
		//we are already looking at new mov flags
		if(!summ_to_previous)
		  {
		    ipath++;
		    su3_put_to_zero(paths[ipath]);
		  }
	      }
	    else su3_put_to_zero(paths[ipath]);
	    
	    //reset the path
	    su3_put_to_id(cur_path);
	  }
	
	//look whether we need to use local or non-local buffer
	su3 *links;
	if(nonloc) links=nonloc_links;
	else       links=(su3*)conf;
	
	//multiply for the link or the link daggered
	if(herm) safe_su3_prod_su3_dag(cur_path,cur_path,links[ilink]);
	else     safe_su3_prod_su3    (cur_path,cur_path,links[ilink]);
      }
    
    //summ last path
    su3_summassign(paths[npaths-1],cur_path);
    
    if(ipath!=npaths-1) crash("finished the movements at path %d while expecting %d",ipath,npaths);
    
    nissa_free(nonloc_links);
  }
  
  void paths_calculation_structure::compute_eo(su3 *paths,quad_su3 **conf)
  {
    //buffer for non local links gathering
    su3 *nonloc_links=nissa_malloc("links",nind_nonloc_links,su3);
    
    //gather non local links
    gather_nonloc_eo(nonloc_links,conf);
    
    //running path
    su3 cur_path;
    su3_put_to_zero(cur_path);
    
    //compute the paths one by one
    int ipath=0;
    for(int imov=0;imov<ntot_mov;imov++)
      {
	int ilink=link_for_movements[imov]>>nposs_path_flags;
	
	//take the flags
	int tag=link_for_movements[imov]%(1<<nposs_path_flags);
	int start=tag&START_PATH_FLAG;
	int herm=tag&DAG_LINK_FLAG;
	int nonloc=tag&NONLOC_LINK_FLAG;
	int summ_to_previous=tag&SUMM_TO_PREVIOUS_PATH_FLAG;
	
	//check if starting a new path
	if(start==1)
	  {
	    //if not the first mov, start the new path
	    if(imov!=0)
	      {
		//summ the ultimated path in its place
		su3_summassign(paths[ipath],cur_path);
		
		//we are already looking at new mov flags
		if(!summ_to_previous)
		  {
		    ipath++;
		    su3_put_to_zero(paths[ipath]);
		  }
	      }
	    else su3_put_to_zero(paths[ipath]);
	    
	    //reset the path
	    su3_put_to_id(cur_path);
	  }
	
	//look whether we need to use local or non-local buffer
	su3 *link;
	if(nonloc) link=nonloc_links+ilink;
	else
	  {
	    int lx=ilink/4,mu=ilink%4;
	    int p=loclx_parity[lx],eo=loceo_of_loclx[lx];
	    
	    link=conf[p][eo]+mu;
	  }
	
	//multiply for the link or the link daggered
	if(herm) safe_su3_prod_su3_dag(cur_path,cur_path,*link);
	else     safe_su3_prod_su3    (cur_path,cur_path,*link);
      }
    
    //summ last path
    su3_summassign(paths[npaths-1],cur_path);
    
    if(ipath!=npaths-1) crash("finished the movements at path %d while expecting %d",ipath,npaths);
    
    nissa_free(nonloc_links);
  }
  
  /////////////////////////////////////////// A SIMPLER APPROACH /////////////////////////////////////////
  
  //initialise to identity
  void init_su3_path(path_drawing_t* c,su3* out)
  {
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      su3_put_to_id(out[ivol.nastyConvert()]);
    NISSA_PARALLEL_LOOP_END;
    coords_t t;
    c->push_back(t);
    
    set_borders_invalid(out);
  }
  
  //compare start and end
  void crash_if_end_diff_from_start(path_drawing_t *c)
  {
    if(c->back()!=c->front())
      crash("(end point,[%d %d %d %d])!=(start_point,[%d,%d,%d,%d])",
	    c->back()[0],c->back()[1],c->back()[2],c->back()[3],
	    c->front()[0],c->front()[1],c->front()[2],c->front()[3]);
  }
  
  //elong backward
  void elong_su3_path_BW(path_drawing_t* c,su3* out,quad_su3* conf,int mu,bool both_sides)
  {
    
    if(both_sides) crash_if_end_diff_from_start(c);
    
    su3_vec_single_shift(out,mu,-1);
    
    if(both_sides)
      {
	NISSA_PARALLEL_LOOP(ivol,0,locVol)
	  {
	    su3 temp;
	    unsafe_su3_prod_su3_dag(temp,out[ivol.nastyConvert()],conf[ivol.nastyConvert()][mu]);
	    unsafe_su3_prod_su3(out[ivol.nastyConvert()],conf[ivol.nastyConvert()][mu],temp);
	  }
	NISSA_PARALLEL_LOOP_END;
      }
    else
      {
	NISSA_PARALLEL_LOOP(ivol,0,locVol)
	  safe_su3_prod_su3_dag(out[ivol.nastyConvert()],out[ivol.nastyConvert()],conf[ivol.nastyConvert()][mu]);
	NISSA_PARALLEL_LOOP_END;
      }
    
    coords_t t(c->back());
    t[mu]--;
    c->push_back(t);
    if(both_sides) c->push_front(t);
    
    THREAD_BARRIER();
  }
  
  //elong forward
  void elong_su3_path_FW(path_drawing_t* c,su3* out,quad_su3* conf,int mu,bool both_sides)
  {
    
    if(both_sides) crash_if_end_diff_from_start(c);
    
    if(both_sides)
      {
	NISSA_PARALLEL_LOOP(ivol,0,locVol)
	  {
	    su3 temp;
	    unsafe_su3_prod_su3(temp,out[ivol.nastyConvert()],conf[ivol.nastyConvert()][mu]);
	    unsafe_su3_dag_prod_su3(out[ivol.nastyConvert()],conf[ivol.nastyConvert()][mu],temp);
	  }
	NISSA_PARALLEL_LOOP_END;
      }
    else
      {
	NISSA_PARALLEL_LOOP(ivol,0,locVol)
	  safe_su3_prod_su3(out[ivol.nastyConvert()],out[ivol.nastyConvert()],conf[ivol.nastyConvert()][mu]);
	NISSA_PARALLEL_LOOP_END;
      }
    THREAD_BARRIER();
    
    coords_t t(c->back());
    t[mu]++;
    c->push_back(t);
    if(both_sides) c->push_front(t);
    
    su3_vec_single_shift(out,mu,+1);
  }
  
  //elong of a certain numer of steps in a certain oriented direction: -1=BW, +1=FW
  void elong_su3_path(path_drawing_t* c,su3* out,quad_su3* conf,int mu,int len,bool both_sides)
  {
    //pointer to avoid branch
    void (*fun)(path_drawing_t*,su3*,quad_su3*,int,bool)=((len<0)?elong_su3_path_BW:elong_su3_path_FW);
    
    //call the appropriate number of times
    for(int l=0;l<abs(len);l++) fun(c,out,conf,mu,both_sides);
  }
  
  //elong a path following a number of macro-steps
  //each step is a pairs consisting of a direction, and length with an orientation
  void elong_su3_path(path_drawing_t* c,su3* out,quad_su3* conf,path_list_steps_t steps)
  {
    for(int i=0;i<(int)steps.size();i++) elong_su3_path(c,out,conf,steps[i]);
  }
}
