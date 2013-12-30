#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "all_to_all.hpp"

namespace nissa
{  
  all_to_all_comm_t::~all_to_all_comm_t()
  {
    nissa_free(list_ranks_to);
    nissa_free(list_ranks_fr);
    nissa_free(in_buf_dest);
    nissa_free(out_buf_source);
    nissa_free(nper_rank_fr);
    nissa_free(nper_rank_to);
    nissa_free(out_buf_off_per_rank);
    nissa_free(in_buf_off_per_rank);
  }

  //find (communicating) the complementary info
  void all_to_all_comm_t::setup_nper_rank_ot_temp(int *nper_rank_ot_temp,int *nper_rank_temp)
  {
    GET_THREAD_ID();
    
    //send and receive the number of elements to communicate with
    if(IS_MASTER_THREAD)
      for(int delta_rank=0;delta_rank<nranks;delta_rank++)
	{
	  int dest_rank=(rank+nranks+delta_rank)%nranks;
	  int recv_rank=(rank+nranks-delta_rank)%nranks;
	  MPI_Sendrecv(nper_rank_temp+dest_rank,1,MPI_INT,dest_rank,0,
		       nper_rank_ot_temp+recv_rank,1,MPI_INT,recv_rank,0,
		       MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
    THREAD_BARRIER();
  }

  //common part of initialization
  void all_to_all_comm_t::common_setup_part1(temp_build_t &build)
  {
    //count the number of different ranks to send and receive from
    int nranks_to_loc=0,nranks_fr_loc=0; //local: structure could be global!
    for(int irank=0;irank<nranks;irank++)
      {
	if(build.nper_rank_to_temp[irank]!=0) nranks_to_loc++;
	if(build.nper_rank_fr_temp[irank]!=0) nranks_fr_loc++;
      }
    nranks_to=nranks_to_loc;
    nranks_fr=nranks_fr_loc;
    list_ranks_to=nissa_malloc("list_ranks_to",nranks_to,int);
    list_ranks_fr=nissa_malloc("list_ranks_fr",nranks_fr,int);

    //store the true list of ranks to communicate with and the number of elements to send to each rank
    nper_rank_to=nissa_malloc("nper_rank_to",nranks_to,int);
    nper_rank_fr=nissa_malloc("nper_rank_fr",nranks_fr,int);
    int irank_to=0,irank_fr=0;
    for(int irank=0;irank<nranks;irank++)
      {
	if(build.nper_rank_to_temp[irank]!=0)
	  {
	    list_ranks_to[irank_to]=irank;
	    nper_rank_to[irank_to]=build.nper_rank_to_temp[irank];
	    build.rank_to_map_list_ranks_to[irank]=irank_to;
	    
	    irank_to++;
	  }
	if(build.nper_rank_fr_temp[irank]!=0)
	  {
	    list_ranks_fr[irank_fr]=irank;
	    nper_rank_fr[irank_fr]=build.nper_rank_fr_temp[irank];
	    build.rank_fr_map_list_ranks_fr[irank]=irank_fr;
	    
	    irank_fr++;
	  }
      }
    
    //count the number of element to receive
    int nel_in_loc=0,nel_out_loc=0;
    for(int irank=0;irank<nranks;irank++)
      {
	nel_in_loc+=build.nper_rank_fr_temp[irank];
	nel_out_loc+=build.nper_rank_to_temp[irank];
      }
    nel_in=nel_in_loc;
    nel_out=nel_out_loc;
    
    //compute the offset and the current indexing element (temporary info)
    out_buf_off_per_rank=nissa_malloc("out_buf_off",nranks_to,int);
    in_buf_off_per_rank=nissa_malloc(" in_buf_off",nranks_fr,int);
    build.out_buf_cur_per_rank=nissa_malloc("out_buf_cur",nranks_to,int);
    build.in_buf_cur_per_rank=nissa_malloc(" in_buf_cur",nranks_fr,int);
    out_buf_off_per_rank[0]=build.out_buf_cur_per_rank[0]=0;
    for(int ilist_rank_to=1;ilist_rank_to<nranks_to;ilist_rank_to++)
      out_buf_off_per_rank[ilist_rank_to]=build.out_buf_cur_per_rank[ilist_rank_to]=
	out_buf_off_per_rank[ilist_rank_to-1]+nper_rank_to[ilist_rank_to-1];
    in_buf_off_per_rank[0]=build.in_buf_cur_per_rank[0]=0;
    for(int ilist_rank_fr=1;ilist_rank_fr<nranks_fr;ilist_rank_fr++)
      in_buf_off_per_rank[ilist_rank_fr]=build.in_buf_cur_per_rank[ilist_rank_fr]=
	in_buf_off_per_rank[ilist_rank_fr-1]+nper_rank_fr[ilist_rank_fr-1];
  }
  
  //explain to each rank how to sort out or fill data from the buff
  //"expl" target the ranks to instruct
  //"note" means the ranks from which learning
  void all_to_all_comm_t::common_setup_part2(int nel_note,
    int *&buf_note,int nranks_note,int *list_ranks_note,int *buf_note_off_per_rank,int *nper_rank_note, 
    int *buf_expl,int nranks_expl,int *list_ranks_expl,int *buf_expl_off_per_rank,int *nper_rank_expl)
  {
    GET_THREAD_ID();
    
    buf_note=nissa_malloc("buf_note",nel_note,int);
    if(IS_MASTER_THREAD)
      {
	MPI_Request req_list[nranks_note+nranks_expl];
	int ireq=0;
	for(int irank_expl=0;irank_expl<nranks_expl;irank_expl++)
	  MPI_Isend(buf_expl+buf_expl_off_per_rank[irank_expl],nper_rank_expl[irank_expl],MPI_INT,
		    list_ranks_expl[irank_expl],909+rank*nranks+list_ranks_expl[irank_expl],cart_comm,&req_list[ireq++]);
	for(int irank_note=0;irank_note<nranks_note;irank_note++)
	  MPI_Irecv(buf_note+buf_note_off_per_rank[irank_note],nper_rank_note[irank_note],MPI_INT,
		    list_ranks_note[irank_note],909+list_ranks_note[irank_note]*nranks+rank,cart_comm,&req_list[ireq++]);
	if(ireq!=nranks_note+nranks_expl) crash("expected %d request, obtained %d",nranks_note+nranks_expl,ireq);
	MPI_Waitall(ireq,req_list,MPI_STATUS_IGNORE);
      }
    THREAD_BARRIER();
  
    //check
    /*
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
    */
  }
  
  //build knowing where to send
  all_to_all_comm_t::all_to_all_comm_t(all_to_all_scattering_list_t &sl)
  {
    setup_knowing_where_to_send(sl);
  }
  void all_to_all_comm_t::setup_knowing_where_to_send(all_to_all_scattering_list_t &sl)
  {
    GET_THREAD_ID();

    temp_build_t build;
    nel_out=sl.size();
    
    //count how many elements to send to each rank
    vector_reset(build.nper_rank_to_temp);
    if(IS_MASTER_THREAD)
      for(int iel_out=0;iel_out<nel_out;iel_out++)
	{
	  int rank_to=sl[iel_out].rank_iel_to.rank;
	  if(rank_to>=nranks||rank_to<0) crash("destination rank %d does not exist!",rank_to);
	  build.nper_rank_to_temp[rank_to]++;
	}
    THREAD_BARRIER();
    
    //count how many elements to receive to each rank
    setup_nper_rank_ot_temp(build.nper_rank_fr_temp,build.nper_rank_to_temp);
    common_setup_part1(build);
    out_buf_source=nissa_malloc("out_buf_source",nel_out,int);

    //save where to store and fill where to load each element in the temporary buffer
    int *in_buf_dest_expl=nissa_malloc("in_buf_dest_expl",nel_out,int);
    if(IS_MASTER_THREAD)
      for(int iel_out=0;iel_out<nel_out;iel_out++)
	{
          int iel_to=sl[iel_out].rank_iel_to.iel,rank_to=sl[iel_out].rank_iel_to.rank;
	  int ilist_rank_to=build.rank_to_map_list_ranks_to[rank_to];
	  int ipos=build.out_buf_cur_per_rank[ilist_rank_to]++;
	  out_buf_source[ipos]=sl[iel_out].iel_fr;
	  in_buf_dest_expl[ipos]=iel_to;
	}
    
    common_setup_part2(nel_in,in_buf_dest,nranks_fr,list_ranks_fr,in_buf_off_per_rank,nper_rank_fr, 
		       in_buf_dest_expl,nranks_to,list_ranks_to,out_buf_off_per_rank,nper_rank_to);
    nissa_free(in_buf_dest_expl);
  }

  //build knowing where to send
  all_to_all_comm_t::all_to_all_comm_t(all_to_all_gathering_list_t &gl)
  {
    setup_knowing_what_to_ask(gl);
  }
  void all_to_all_comm_t::setup_knowing_what_to_ask(all_to_all_gathering_list_t &gl)
  {
    GET_THREAD_ID();

    temp_build_t build;
    nel_in=gl.size();
        
    //count how many elements to send to each rank
    vector_reset(build.nper_rank_fr_temp);
    if(IS_MASTER_THREAD)
      for(all_to_all_gathering_list_t::iterator it=gl.begin();it!=gl.end();it++)
	{
	  int rank_fr=it->rank_iel_fr.rank;
	  if(rank_fr>=nranks||rank_fr<0) crash("source rank %d does not exist!",rank_fr);
	  build.nper_rank_fr_temp[rank_fr]++;
	}
    THREAD_BARRIER();
    
    //count how many elements to receive from each rank
    setup_nper_rank_ot_temp(build.nper_rank_to_temp,build.nper_rank_fr_temp);
    common_setup_part1(build);
    in_buf_dest=nissa_malloc("in_buf_dest",nel_in,int);

    //save the explenations to each rank on how to fill outbuffers
    int *out_buf_source_expl=nissa_malloc("out_buf_source_expl",nel_in,int);
    if(IS_MASTER_THREAD)
      for(all_to_all_gathering_list_t::iterator it=gl.begin();it!=gl.end();it++)
	{
	  int iel_fr=it->rank_iel_fr.iel,rank_fr=it->rank_iel_fr.rank;
	  int ilist_rank_fr=build.rank_fr_map_list_ranks_fr[rank_fr];
	  int ipos=build.in_buf_cur_per_rank[ilist_rank_fr]++;
	  in_buf_dest[ipos]=it->iel_to;
	  out_buf_source_expl[ipos]=iel_fr;
	}
    
    common_setup_part2(nel_out,out_buf_source,nranks_to,list_ranks_to,out_buf_off_per_rank,nper_rank_to, 
		       out_buf_source_expl,nranks_fr,list_ranks_fr,in_buf_off_per_rank,nper_rank_fr);
    nissa_free(out_buf_source_expl);
  }
  
  //perform the remapping
  void all_to_all_comm_t::communicate(void *out,void *in,int bps,void *ext_out_buf,void *ext_in_buf)
  {
    GET_THREAD_ID();
    
    //allocate a buffer where to repack data
    char *out_buf=(ext_out_buf==NULL)?nissa_malloc("out_buf",nel_out*bps,char):(char*)ext_out_buf;
    char *in_buf=(ext_in_buf==NULL)?nissa_malloc("in_buf",nel_in*bps,char):(char*)ext_in_buf;
    
    //copy data on the out-going buffer
    NISSA_PARALLEL_LOOP(iel_out,0,nel_out)
      memcpy(out_buf+iel_out*bps,(char*)in+out_buf_source[iel_out]*bps,bps);
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

    if(ext_out_buf==NULL) nissa_free(out_buf);
    if(ext_in_buf==NULL) nissa_free(in_buf);
  }
}
