#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <string.h>

#include "../base/global_variables.h"
#include "../base/vectors.h"
#include "../base/debug.h"
#include "../geometry/geometry_lx.h"

//remap a vector across all the ranks
void remap_vector(char *out,char *in,coords *xto,coords *xfr,int bps)
{
  //allocate a buffer where to repack data
  char *out_buf=nissa_malloc("out_buf",loc_vol*(bps+sizeof(int)),char);
  char *in_buf=nissa_malloc("in_buf",loc_vol*(bps+sizeof(int)),char);
  
  //allocate addresses, counting and reordering vectors
  int *rank_fr=nissa_malloc("rank_from",loc_vol,int);
  int *rank_to=nissa_malloc("rank_to",loc_vol,int);
  int *ivol_to=nissa_malloc("ivol_to",loc_vol,int);
  int *nfr=nissa_malloc("nfr",nissa_nranks,int);
  int *nto=nissa_malloc("nto",nissa_nranks,int);
  
  //reset the count of data from each node
  memset(nfr,0,sizeof(int)*nissa_nranks);
  memset(nto,0,sizeof(int)*nissa_nranks);
  
  //scan all local sites to see where to send and from where to expect data
  nissa_loc_vol_loop(ivol)
    {
      //compute destination site and destination rank
      get_loclx_and_rank_of_coord(&(ivol_to[ivol]),&(rank_to[ivol]),xto[ivol]);
      if(rank_to[ivol]>=nissa_nranks||rank_to[ivol]<0)
	crash("index %d should go to exceeding rank %d",ivol,rank_to[ivol]);
      
      //count data to be sent to each rank
      nto[rank_to[ivol]]++;
    }
  
  //if we were provided with info on what to expect
  if(xfr!=NULL)
    nissa_loc_vol_loop(ivol)
    {
      //compute origin rank
      rank_fr[ivol]=rank_hosting_site_of_coord(xfr[ivol]);
      //count data to be received from each rank
      nfr[rank_fr[ivol]]++;
    }
  else
    for(int delta_rank=0;delta_rank<nissa_nranks;delta_rank++)
      {
	//we must sort it out: send to +delta, receive from -delta
	int dest_rank=(rank+nissa_nranks+delta_rank)%nissa_nranks;
	int recv_rank=(rank+nissa_nranks-delta_rank)%nissa_nranks;
	MPI_Sendrecv(nto+dest_rank,1,MPI_INT,dest_rank,0,
		     nfr+recv_rank,1,MPI_INT,recv_rank,0,
		     MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
  
  //allocate starting position of out going buffer
  char **out_buf_rank=nissa_malloc("out_buf_rank",nissa_nranks,char*);
  out_buf_rank[0]=out_buf;
  for(int irank=1;irank<nissa_nranks;irank++)
    out_buf_rank[irank]=out_buf_rank[irank-1]+nto[irank-1]*(bps+sizeof(int));
  
  //copy data on the out-going buffer
  nissa_loc_vol_loop(ivol)
    {
      int trank=rank_to[ivol];
      memcpy(out_buf_rank[trank],in+ivol*bps,bps);
      out_buf_rank[trank]+=bps;
      (*((int*)(out_buf_rank[trank])))=ivol_to[ivol];
      out_buf_rank[trank]+=sizeof(int);
    }
  
  //now open communication from each rank which need to send data
  //and to each rank which needs to receive data
  MPI_Request req_list[2*nissa_nranks];
  int ireq=0;
  char *send_buf=out_buf; //init send pointer
  char *recv_buf=in_buf;  //init recv pointer
  char *internal_send_buf=NULL;
  char *internal_recv_buf=NULL;
  for(int irank=0;irank<nissa_nranks;irank++)
    {
      if(rank!=irank)
	{
	  //start communications with other ranks
	  if(nfr[irank]!=0)
	    MPI_Irecv(recv_buf,nfr[irank]*(bps+sizeof(int)),MPI_CHAR,irank,909+irank*nissa_nranks+rank,cart_comm,&req_list[ireq++]);
	  if(nto[irank]!=0)
	    MPI_Isend(send_buf,nto[irank]*(bps+sizeof(int)),MPI_CHAR,irank,909+rank*nissa_nranks+irank,cart_comm,&req_list[ireq++]);
	}
      else
	{
	  //save pointer to internal communications
	  internal_send_buf=send_buf;
	  internal_recv_buf=recv_buf;
	}
      recv_buf+=nfr[irank]*(bps+sizeof(int));
      send_buf+=nto[irank]*(bps+sizeof(int));
    }
  if(recv_buf!= in_buf+(sizeof(int)+bps)*loc_vol) crash("Mismatch in receiving buffer");
  if(send_buf!=out_buf+(sizeof(int)+bps)*loc_vol) crash("Mismatch in sending buffer");
  
  //now, while waiting for all requests to be accomplished, copy internal data
  memcpy(internal_recv_buf,internal_send_buf,(sizeof(int)+bps)*nfr[rank]);
  
  //wait to accomplish all the requests
  MPI_Waitall(ireq,req_list,MPI_STATUS_IGNORE);
  
  //now sort out data from the incoming buffer
  nissa_loc_vol_loop(ivol)
    {
      char *start_in_buf=in_buf+ivol*(bps+sizeof(int));
      int to=(*((int*)(start_in_buf+bps)));
      memcpy(out+to*bps,start_in_buf,bps);
    }
  
  //invalidate borders
  set_borders_invalid(out);
  
  //free vectors
  
  nissa_free(out_buf_rank);
  
  nissa_free(nfr);
  nissa_free(nto);
  
  nissa_free(ivol_to);
  nissa_free(rank_to);
  nissa_free(rank_fr);
  
  nissa_free(out_buf);
  nissa_free(in_buf);
}

