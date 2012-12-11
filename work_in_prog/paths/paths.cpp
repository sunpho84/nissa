#include "nissa.h"
#include <mpi.h>

int T=32,L=16;
int debug=0;

struct ind_links_vector
{
  //int nlinks
  int buf_size;
  int nind_links;
  int *ind_links;
  char *links_buf;
  ind_links_vector() : nind_links(0),buf_size(0),ind_links(NULL) {
    links_buf=nissa_malloc("links_buf",glb_vol*4/8+1,char);
    vector_reset(links_buf);
  }
  int search(int t);
  int search(int in,int mu) {return search(combine(in,mu));}
  int append(int t);
  int append(int in,int mu) {return append(combine(in,mu));}
  void finish() {if(links_buf==NULL) nissa_free(links_buf);}
  ~ind_links_vector() {if(links_buf!=NULL) nissa_free(links_buf);if(ind_links!=NULL) free((void*)ind_links);}
  int combine(int in,int mu) {if(in<0||in>=glb_vol) crash("x=%d must be in [0,%d]",in,glb_vol-1);if(mu<0||mu>=4) crash("mu=%d, must be in [0,3]");return in*4+mu;}
};

int ind_links_vector::search(int t)
{
  int i=0;
  while(i<nind_links && ind_links[i]!=t) i++;
  
  return i;
}

int ind_links_vector::append(int t)
{
  if(debug) master_printf("Appending %d\n",t);

  if(links_buf[t/8]&(1<<(t%8)))
    {
      int pos=search(t);
      if(pos==nind_links) crash("something went wrong");
      return pos;
    }
  else
    {
      if(nind_links>=buf_size)
	{
	  buf_size=2*buf_size+1;
	  ind_links=(int*)realloc((void*)ind_links,sizeof(int)*buf_size);
	  printf("Reallocating: %d\n",buf_size);
	}
      ind_links[nind_links]=t;
      nind_links++;
      links_buf[t/8]|=(1<<(t%8));

      return nind_links-1;
    }
}

//NB: the list of movement is larger than the list of indipendent links!
struct paths_calculation_structure
{
  //initialization
  paths_calculation_structure(int npaths,int ntot_mov) : npaths(npaths),ntot_mov(ntot_mov),cur_path(0),cur_mov(0) {
    ind_links_list=new ind_links_vector;
    master_printf("initializing a new path calculation structure with %d movements and %d paths\n",ntot_mov,npaths);
    link_for_movements=nissa_malloc("link_for_movements",ntot_mov,int);
    //prepare the list of communication structures
    nranks_to_send=nranks_to_recv=0;
    ranks_to_send=ranks_to_recv=NULL;
    masks_to_recv=masks_to_send=NULL;
  };
  ~paths_calculation_structure() {
    delete ind_links_list;
    nissa_free(link_for_movements);
    
    for(int irank=0;irank<nranks_to_send;irank++) MPI_Type_free(masks_to_send+irank);
    for(int irank=0;irank<nranks_to_recv;irank++) MPI_Type_free(masks_to_recv+irank);
    
    free((void*)ranks_to_send);
    free((void*)ranks_to_recv);
    free((void*)masks_to_send);
    free((void*)masks_to_recv);
  }
  
  //parameters defining the set of paths
  int npaths,ntot_mov;
  int *link_for_movements;
  ind_links_vector *ind_links_list;
  
  //current global movement (link), path and position
  int cur_path,cur_mov;
  int pos;
  
  //relevant for MPI part
  int nranks_to_send,*ranks_to_send;
  int nranks_to_recv,*ranks_to_recv;
  MPI_Datatype *masks_to_recv;
  MPI_Datatype *masks_to_send;
  
  //commands
  void start_new_path_from_loclx(int lx) {
    pos=glblx_of_loclx[lx];
    link_for_movements[cur_mov]=1;
    if(debug) master_printf("Starting a new path from local point %d, global point: %d\n",lx,pos);}
  void switch_to_next_step() {
    cur_mov++;
    if(cur_mov>ntot_mov) crash("exceded the number of allocatec movements, %d",ntot_mov);
    link_for_movements[cur_mov]=0;}
  void move_forward(int mu);
  void move_backward(int mu);
  void stop_current_path() {
    cur_path++;
    if(cur_path>npaths) crash("exceded the number of allocated paths, %d",npaths);
  }
  void finished_last_path();
  void communicate(su3 *paths,quad_su3 *conf);
  void compute(su3 *paths,quad_su3 *conf);
 
private:
  //avoid bare initialization without specification of nel
  paths_calculation_structure();
  //setu the communication buffers
  void setup_sender_receivers();
};

//add a forward move
void paths_calculation_structure::move_forward(int mu)
{
  //check not to have passed the max number of steps
  if(cur_mov==ntot_mov) crash("exceded the number of allocated movements, %d",ntot_mov);
  //find global pos
  int n=glblx_neighup(pos,mu);
  //shift by 2 not to overwrite tags
  int ul=ind_links_list->append(pos,mu);
  link_for_movements[cur_mov]+=(ul<<2);
  //set new pos
  if(debug) master_printf("Moved forward from %d in the direction %d to %d, mov: %d, link: %d, tag: %d\n",pos,mu,n,cur_mov,ul,link_for_movements[cur_mov]%4);
  pos=n;
  //switch to next step
  switch_to_next_step();
}

//add a backward move
void paths_calculation_structure::move_backward(int mu)
{
  //check not to have passed the max number of steps
  if(cur_mov==ntot_mov) crash("exceded the number of allocated movements, %d",ntot_mov);
  //find global pos
  int n=glblx_neighdw(pos,mu);
  //shift by 2 not to overwrite tags, and append the "dagger" flag
  int ul=ind_links_list->append(n,mu);
  link_for_movements[cur_mov]+=(ul<<2)+2;
  //set new pos
  if(debug) master_printf("Moved backward from %d in the direction %d to %d, mov: %d, link: %d, tag: %d\n",pos,mu,n,cur_mov,ul,link_for_movements[cur_mov]%4);
  pos=n;
  //switch to next step
  switch_to_next_step();
}

//finish settings the paths, setup the send and receiver
void paths_calculation_structure::finished_last_path()
{
  if(cur_path!=npaths) crash("finished the path list at path %d while it was initialized for %d",cur_path,npaths);
  if(cur_mov!=ntot_mov) crash("finished the path list at mov %d while it was initialized for %d",cur_mov,ntot_mov);
  
  master_printf("The calculation of all the paths involves %d indipendent links\n",ind_links_list->nind_links);
  setup_sender_receivers();
}

//setup the sender and receiver buffers, finding which ranks are involved
void paths_calculation_structure::setup_sender_receivers()
{
  //loop over the ranks displacement, counting local
  for(int delta_rank=1;delta_rank<nissa_nranks;delta_rank++)
    {
      int send_to_rank=(rank+delta_rank)%nissa_nranks;
      int recv_fr_rank=(rank+nissa_nranks-delta_rank)%nissa_nranks;
      
      //counts the number of links to receive
      int nlinks_to_recv=0;
      for(int i=0;i<ind_links_list->nind_links;i++)
	if(rank_hosting_glblx(ind_links_list->ind_links[i]>>2)==recv_fr_rank)
	  nlinks_to_recv++;
      
      //send this piece of info and receive the number of links to send
      int nlinks_to_send=0;
      MPI_Sendrecv((void*)&(nlinks_to_recv),1,MPI_INT,recv_fr_rank,recv_fr_rank*nissa_nranks+rank,
		   (void*)&(nlinks_to_send),1,MPI_INT,send_to_rank,rank*nissa_nranks+send_to_rank,
		   MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      
      if(debug)
	{
	  //print the number of links to ask to and asked by other ranks
	  printf("rank %d, recv from rank %d: %d links, send to rank %d: %d links\n",rank,recv_fr_rank,nlinks_to_recv,send_to_rank,nlinks_to_send);
	}
      
      //allocate a buffer where to store the list of links to receive, and to ask (differs for ordering)
      int *links_to_ask=nissa_malloc("links_to_ask",nlinks_to_recv,int);
      int *links_to_recv=nissa_malloc("links_to_recv",nlinks_to_recv,int);
      int ilink=0;
      for(int i=0;i<ind_links_list->nind_links;i++)
	{
	  int t=ind_links_list->ind_links[i];
	  int gx=t>>2;
	  int mu=t%4;
	  
	  if(debug)
	    {
	      int flag=t%4;
	      if(delta_rank==1) printf("rank %d, link %d, flag %d\n",rank,i,flag);
	    }
	  
	  //get lx and rank hosting the site
	  int lx,rx;
	  get_loclx_and_rank_of_glblx(&lx,&rx,gx);
	  
	  //copy in the list if appropriate rank
	  if(rx==recv_fr_rank)
	    {
	      links_to_ask[ilink]=(lx<<2)+mu;
	      links_to_recv[ilink]=i;
	      ilink++;
	    }
	}
      if(debug)
	{
	  printf("rank %d, total link to recv: %d\n",rank,ilink);
	  for(int ilink=0;ilink<nlinks_to_recv;ilink++)
	    printf("rank %d will recv from %d as link %d/%d link: %d\n",rank,recv_fr_rank,ilink,nlinks_to_recv,links_to_recv[ilink]);      
	}
      
      //allocate the list of link to send
      int *links_to_send=nissa_malloc("links_to_send",nlinks_to_send,int);
      
      //send this piece of info
      MPI_Sendrecv((void*)links_to_ask, nlinks_to_recv,MPI_INT,recv_fr_rank,recv_fr_rank*nissa_nranks+rank,
                   (void*)links_to_send,nlinks_to_send,MPI_INT,send_to_rank,rank*nissa_nranks+send_to_rank,
                   MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      nissa_free(links_to_ask);
      
      if(debug)
	for(int ilink=0;ilink<nlinks_to_send;ilink++)
	  printf("rank %d will send to %d as link %d/%d link: %d\n",rank,send_to_rank,ilink,nlinks_to_send,links_to_send[ilink]);
      
      //store the sending rank id, make space for it, create the type and commit
      if(nlinks_to_send!=0)
	{
	  //prepare the space
	  nranks_to_send++;
	  ranks_to_send=(int*)realloc(ranks_to_send,sizeof(int)*nranks_to_send);
	  masks_to_send=(MPI_Datatype*)realloc(masks_to_send,sizeof(MPI_Datatype)*nranks_to_send);
	  ranks_to_send[nranks_to_send-1]=send_to_rank;
	  
	  //prepare the repetition list
	  int *single=nissa_malloc("single",nlinks_to_send,int);
	  for(int ilink=0;ilink<nlinks_to_send;ilink++)
	    {
	      single[ilink]=1;
	      if(debug) printf("rank %d sending link: %d\n",rank,links_to_send[ilink]);
	    }
	  
	  //create the datatype
	  MPI_Type_indexed(nlinks_to_send,single,links_to_send,MPI_SU3,masks_to_send+nranks_to_send-1);
	  MPI_Type_commit(masks_to_send+nranks_to_send-1);
	  
	  if(debug) printf("rank %d, mask to send %d created, nlink: %d\n",rank,nranks_to_send-1,nlinks_to_send);
	}
      
      //store the receiving rank id, make space for it, create the type and commit
      if(nlinks_to_recv!=0)
	{
	  //prepare the space
	  nranks_to_recv++;
	  ranks_to_recv=(int*)realloc(ranks_to_recv,sizeof(int)*nranks_to_recv);
	  masks_to_recv=(MPI_Datatype*)realloc(masks_to_recv,sizeof(MPI_Datatype)*nranks_to_recv);
	  ranks_to_recv[nranks_to_recv-1]=recv_fr_rank;
	  
	  //prepare the repetition list
          int *single=nissa_malloc("single",nlinks_to_recv,int);
          for(int ilink=0;ilink<nlinks_to_recv;ilink++)
	    {
	      single[ilink]=1;
	      if(debug) printf("rank %d receiving link: %d\n",rank,links_to_recv[ilink]);
	    }
	  
	  //create the datatype
	  MPI_Type_indexed(nlinks_to_recv,single,links_to_recv,MPI_SU3,masks_to_recv+nranks_to_recv-1);
	  MPI_Type_commit(masks_to_recv+nranks_to_recv-1);
	  
	  nissa_free(single);
	}
      nissa_free(links_to_recv);
      nissa_free(links_to_send);
    }
}

//collect all the indipendent links entering into the calculations into "links"
void paths_calculation_structure::communicate(su3 *links,quad_su3 *conf)
{
   //list of request
  MPI_Request request[nranks_to_send+nranks_to_recv];
  int irequest=0;
  
  //open sending communications
  for(int isend=0;isend<nranks_to_send;isend++)
    {
      MPI_Isend((void*)conf,1,masks_to_send[isend],ranks_to_send[isend],ranks_to_send[isend]*nissa_nranks+rank,cart_comm,request+irequest++);
      int size;
      MPI_Type_size(masks_to_send[isend],&size);
      if(debug) printf("Rank %d sending to %d, mask %d, size: %d\n",rank,ranks_to_send[isend],isend,size);
    }
  
  //open receiving communications
  for(int irecv=0;irecv<nranks_to_recv;irecv++)
    {  
      MPI_Irecv((void*)links,1,masks_to_recv[irecv],ranks_to_recv[irecv],rank*nissa_nranks+ranks_to_recv[irecv],cart_comm,request+irequest++);
      int size;
      MPI_Type_size(masks_to_recv[irecv],&size);
      if(debug) printf("Rank %d receiving from %d, mask %d, size: %d\n",rank,ranks_to_recv[irecv],irecv,size);
    }
  
  //in the while, copy internal links
  for(int i=0;i<ind_links_list->nind_links;i++)
    {
      int t=ind_links_list->ind_links[i];
      int gx=t>>2;
      int mu=t%4;
      
      //get lx and rank hosting the site
      int lx,rx;
      get_loclx_and_rank_of_glblx(&lx,&rx,gx);
      
      //copy in the list if appropriate rank
      if(rx==rank) su3_copy(links[i],conf[lx][mu]);
    }
  
  //wait communications to finish
  MPI_Waitall(nranks_to_send+nranks_to_recv,request,MPI_STATUS_IGNORE);
  
  if(debug)
    {
      int ind[2]={1,1024};
      if(rank==1)
	for(int i=0;i<2;i++)
	  {
	    printf("Link to send: %d\n",i);
	    su3_print(((su3*)conf)[ind[i]]);
	  }
      fflush(stdout);
      sleep(2);
      MPI_Barrier(MPI_COMM_WORLD);
      
      if(rank==0)
	for(int i=0;i<3;i++)
	  {
	    printf("Link received: %d\n",i);
	    su3_print(links[i]);
	  }
      fflush(stdout);
      sleep(2);
      MPI_Barrier(MPI_COMM_WORLD);
    }
}

void paths_calculation_structure::compute(su3 *paths,quad_su3 *conf)
{
  //buffer
  su3 *links=nissa_malloc("links",ind_links_list->nind_links,su3);
  
  //communicate
  communicate(links,conf);
  
  //compute the paths one by one
  int ipath=0;
  for(int imov=0;imov<ntot_mov;imov++)
    {
      int ilink=link_for_movements[imov]>>2;
      int tag=link_for_movements[imov]%4;
      int start=tag%2;
      int herm=tag-start;
      if(debug) printf("mov %d, ind link: %d, tag: %d, herm: %d, start: %d\n",imov,ilink,tag,herm,start);
      
      if(start==1)
	{
	  //if not the first mov, start the new path
	  if(imov!=0) ipath++;
	  su3_put_to_id(paths[ipath]);
	}
      
      //multiply for the link or the link daggered
      if(herm) safe_su3_prod_su3_dag(paths[ipath],paths[ipath],links[ilink]);
      else     safe_su3_prod_su3    (paths[ipath],paths[ipath],links[ilink]);
    }
  
  nissa_free(links);
}

int main(int narg,char **arg)
{
  init_nissa();

  init_grid(T,L);

  quad_su3 *conf=nissa_malloc("e",loc_vol,quad_su3);
  read_ildg_gauge_conf(conf,"/Users/francesco/Prace/Confs/conf.0100");

  /////////////////////////
  
  paths_calculation_structure *a=new paths_calculation_structure(loc_vol,loc_vol*4);

  nissa_loc_vol_loop(ivol)
    {
      if(100000*ivol%loc_vol==0) printf("rank %d setting path %d/%d\n",rank,ivol,loc_vol);
      a->start_new_path_from_loclx(ivol);
      a->move_forward(0);
      a->move_forward(1);
      a->move_backward(0);
      a->move_backward(1);
      a->stop_current_path();
    }

  a->finished_last_path();

  su3 *paths=nissa_malloc("paths",loc_vol,su3);
  a->compute(paths,conf);
  
  nissa_free(conf);
  nissa_free(paths);
    
  delete a;
  /*
    link_id_vector a;
    printf("%d\n",a.append(10,2));
    printf("%d\n",a.append(10,1));
    printf("%d\n",a.append(10,2));
    printf("%d\n",a.append(10,1));
  */
  
  /////////////////////////

  close_nissa();
  
  return 0;
}
