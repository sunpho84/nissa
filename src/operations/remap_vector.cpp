#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  
  vector_remap_t::~vector_remap_t()
  {
    nissa_free(list_ranks_to);
    nissa_free(list_ranks_fr);
    nissa_free(in_buf_dest);
    nissa_free(out_buf_pos);
    nissa_free(nper_rank_fr);
    nissa_free(nper_rank_to);
    nissa_free(out_buf_off_per_rank);
    nissa_free(in_buf_off_per_rank);
  }
  
  //constructor
  vector_remap_t::vector_remap_t(int nel_out,void (*index)(int &irank_to,int &iel_to,int iel_fr,void *pars),void *pars) : nel_out(nel_out)
  {
    GET_THREAD_ID();
    
    //count how many elements to send to each rank
    int *nper_rank_to_temp=nissa_malloc("nper_rank_to_temp",nranks,int);
    vector_reset(nper_rank_to_temp);
    if(IS_MASTER_THREAD)
      for(int iel_out=0;iel_out<nel_out;iel_out++)
	{
	  int rank_to,iel_to;
	  index(rank_to,iel_to,iel_out,pars);
	  if(rank_to>=nranks||rank_to<0) crash("destination rank %d does not exist!",rank_to);
	  nper_rank_to_temp[rank_to]++;
	}
    THREAD_BARRIER();
    
    //send and receive the number of elements to receive
    int *nper_rank_fr_temp=nissa_malloc("nper_rank_fr_temp",nranks,int);
    if(IS_MASTER_THREAD)
      for(int delta_rank=0;delta_rank<nranks;delta_rank++)
	{
	  int dest_rank=(rank+nranks+delta_rank)%nranks;
	  int recv_rank=(rank+nranks-delta_rank)%nranks;
	  MPI_Sendrecv(nper_rank_to_temp+dest_rank,1,MPI_INT,dest_rank,0,
		       nper_rank_fr_temp+recv_rank,1,MPI_INT,recv_rank,0,
		       MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
    THREAD_BARRIER();
    
    //count the number of different ranks to send and receive from
    int nranks_to_loc=0,nranks_fr_loc=0; //local: structure could be global!
    for(int irank=0;irank<nranks;irank++)
      {
	if(nper_rank_to_temp[irank]!=0) nranks_to_loc++;
	if(nper_rank_fr_temp[irank]!=0) nranks_fr_loc++;
      }
    nranks_to=nranks_to_loc;
    nranks_fr=nranks_fr_loc;
    list_ranks_to=nissa_malloc("list_ranks_to",nranks_to,int);
    list_ranks_fr=nissa_malloc("list_ranks_fr",nranks_fr,int);

    //store the true list of ranks to communicate with and the number of elements to send to each rank
    nper_rank_to=nissa_malloc("nper_rank_to",nranks_to,int);
    nper_rank_fr=nissa_malloc("nper_rank_fr",nranks_fr,int);
    int irank_to=0,irank_fr=0;
    std::map<int,int> rank_to_map_list_ranks_to,rank_fr_map_list_ranks_fr; //later use
    for(int irank=0;irank<nranks;irank++)
      {
	if(nper_rank_to_temp[irank]!=0)
	  {
	    list_ranks_to[irank_to]=irank;
	    nper_rank_to[irank_to]=nper_rank_to_temp[irank];
	    rank_to_map_list_ranks_to[irank]=irank_to;
	    
	    irank_to++;
	  }
	if(nper_rank_fr_temp[irank]!=0)
	  {
	    list_ranks_fr[irank_fr]=irank;
	    nper_rank_fr[irank_fr]=nper_rank_fr_temp[irank];
	    rank_fr_map_list_ranks_fr[irank]=irank_fr;
	    
	    irank_fr++;
	  }
      }
    
    //count the number of element to receive
    int nel_in_loc=0;
    for(int irank=0;irank<nranks;irank++)
      nel_in_loc+=nper_rank_fr_temp[irank];
    nel_in=nel_in_loc;
    
    //compute the offset and the current indexing element (temporary info)
    out_buf_off_per_rank=nissa_malloc("out_buf_off",nranks_to,int);
    in_buf_off_per_rank=nissa_malloc(" in_buf_off",nranks_fr,int);
    int *out_buf_cur_per_rank=nissa_malloc("out_buf_cur",nranks_to,int);
    int  *in_buf_cur_per_rank=nissa_malloc(" in_buf_cur",nranks_fr,int);
    out_buf_off_per_rank[0]=out_buf_cur_per_rank[0]=0;
    for(int ilist_rank_to=1;ilist_rank_to<nranks_to;ilist_rank_to++)
      out_buf_off_per_rank[ilist_rank_to]=out_buf_cur_per_rank[ilist_rank_to]=
	out_buf_off_per_rank[ilist_rank_to-1]+nper_rank_to[ilist_rank_to-1];
    in_buf_off_per_rank[0]=in_buf_cur_per_rank[0]=0;
    for(int ilist_rank_fr=1;ilist_rank_fr<nranks_fr;ilist_rank_fr++)
      in_buf_off_per_rank[ilist_rank_fr]=in_buf_cur_per_rank[ilist_rank_fr]=
	in_buf_off_per_rank[ilist_rank_fr-1]+nper_rank_fr[ilist_rank_fr-1];
    
    //save where to store and fill where to  load each element in the temporary buffer
    out_buf_pos=nissa_malloc("out_buf_pos",nel_out,int);
    int *in_buf_dest_expl=nissa_malloc("in_buf_dest_expl",nel_out,int);
    if(IS_MASTER_THREAD)
      for(int iel_out=0;iel_out<nel_out;iel_out++)
	{
	  int rank_to,iel_to;
	  index(rank_to,iel_to,iel_out,pars);
	  
	  int ilist_rank_to=rank_to_map_list_ranks_to[rank_to];
	  int ipos=out_buf_cur_per_rank[ilist_rank_to]++;
	  out_buf_pos[iel_out]=ipos;
	  in_buf_dest_expl[ipos]=iel_to;
	}
    
    //explain to each rank how to sort out data from the box
    in_buf_dest=nissa_malloc("in_buf_dest",nel_in,int);
    if(IS_MASTER_THREAD)
      {
	MPI_Request req_list[nranks_to+nranks_fr];
	int ireq=0;
	for(int irank_to=0;irank_to<nranks_to;irank_to++)
	  MPI_Isend(in_buf_dest_expl+out_buf_off_per_rank[irank_to],nper_rank_to[irank_to],MPI_INT,
		    list_ranks_to[irank_to],909+rank*nranks+list_ranks_to[irank_to],cart_comm,&req_list[ireq++]);
	for(int irank_fr=0;irank_fr<nranks_fr;irank_fr++)
	  MPI_Irecv(in_buf_dest+in_buf_off_per_rank[irank_fr],nper_rank_fr[irank_fr],MPI_INT,
		    list_ranks_fr[irank_fr],909+list_ranks_fr[irank_fr]*nranks+rank,cart_comm,&req_list[ireq++]);
	if(ireq!=nranks_to+nranks_fr) crash("expected %d request, obtained %d",nranks_to+nranks_fr,ireq);
	MPI_Waitall(ireq,req_list,MPI_STATUS_IGNORE);
      }
    THREAD_BARRIER();
    
    //check
    int *in_buf_dest_check=nissa_malloc("in_buf_dest_check",nel_in,int);
    vector_reset(in_buf_dest_check);
    if(IS_MASTER_THREAD)
      for(int iel_in=0;iel_in<nel_in;iel_in++)
	{
	  int idest=in_buf_dest[iel_in];
	  if(idest<0||idest>=nel_in) crash("in_buf_dest[%d] point to %d not in the range [0,nel_in=%d)",
					   iel_in,idest,nel_in);
	  if(in_buf_dest_check[idest]++==1) crash("multiple assignement of %d",idest);
	}
    THREAD_BARRIER();
    nissa_free(in_buf_dest_check);
    
    nissa_free(out_buf_cur_per_rank);
    nissa_free(in_buf_cur_per_rank);
    nissa_free(in_buf_dest_expl);
    nissa_free(nper_rank_to_temp);
    nissa_free(nper_rank_fr_temp);
  }
  
  //perform the remapping
  void vector_remap_t::remap(void *out,void *in,int bps)
  {
    GET_THREAD_ID();
    
    //allocate a buffer where to repack data
    char *out_buf=nissa_malloc("out_buf",nel_out*bps,char);
    char *in_buf=nissa_malloc("in_buf",nel_in*bps,char);
    
    //copy data on the out-going buffer
    NISSA_PARALLEL_LOOP(iel_out,0,nel_out)
      memcpy(out_buf+out_buf_pos[iel_out]*bps,(char*)in+iel_out*bps,bps);
    THREAD_BARRIER();
    
    if(IS_MASTER_THREAD)
      {
	MPI_Request req_list[nranks_to+nranks_fr];
	int ireq=0;
	for(int irank_fr=0;irank_fr<nranks_fr;irank_fr++)
	  MPI_Irecv(in_buf+in_buf_off_per_rank[irank_fr]*bps,nper_rank_fr[irank_fr]*bps,MPI_CHAR,
		    list_ranks_fr[irank_fr],909+list_ranks_fr[irank_fr]*nranks+rank,cart_comm,&req_list[ireq++]);
	for(int irank_to=0;irank_to<nranks_to;irank_to++)
	  MPI_Isend(out_buf+out_buf_off_per_rank[irank_to]*bps,nper_rank_to[irank_to]*bps,MPI_CHAR,
		    list_ranks_to[irank_to],909+rank*nranks+list_ranks_to[irank_to],cart_comm,&req_list[ireq++]);
      	if(ireq!=nranks_to+nranks_fr) crash("expected %d request, obtained %d",nranks_to+nranks_fr,ireq);
	MPI_Waitall(ireq,req_list,MPI_STATUS_IGNORE);
      }
    THREAD_BARRIER();
      
    //sort out data from the incoming buffer
    NISSA_PARALLEL_LOOP(iel_in,0,nel_in)
      memcpy((char*)out+in_buf_dest[iel_in]*bps,in_buf+iel_in*bps,bps);
    
    set_borders_invalid(out);

    nissa_free(out_buf);
    nissa_free(in_buf);
  }
  
  void remap_vector(char *out,char *in,coords *xto,coords *xfr,int bps)
  {
  }
}
