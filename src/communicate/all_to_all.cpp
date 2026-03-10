#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <string.h>

#include "base/debug.hpp"
#include "base/memory_manager.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

#include "all_to_all.hpp"

namespace nissa
{
  temp_build_t::temp_build_t()
  {
    nper_rank_to_temp.resize(nranks);
    nper_rank_fr_temp.resize(nranks);
  }
  
  //find (communicating) the complementary info
  void all_to_all_comm_t::setup_nper_rank_other_temp(std::vector<int>& nper_rank_other_temp,
						     const std::vector<int>& nper_rank_temp)
  {
    //send and receive the number of elements to communicate with
    for(int delta_rank=0;delta_rank<nranks;delta_rank++)
      {
	const int dest_rank=(rank+nranks+delta_rank)%nranks;
	const int recv_rank=(rank+nranks-delta_rank)%nranks;
	MPI_Sendrecv(&nper_rank_temp[dest_rank],1,MPI_INT,dest_rank,0,
		     &nper_rank_other_temp[recv_rank],1,MPI_INT,recv_rank,0,
		     MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
    
    //send and receive the number of elements to communicate with
    for(int delta_rank=0;delta_rank<nranks;delta_rank++)
      {
	const int dest_rank=(rank+nranks+delta_rank)%nranks;
	const int recv_rank=(rank+nranks-delta_rank)%nranks;
	MPI_Sendrecv(&nper_rank_temp[dest_rank],1,MPI_INT,dest_rank,0,
		     &nper_rank_other_temp[recv_rank],1,MPI_INT,recv_rank,0,
		     MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
    VERBOSITY_LV3_MASTER_PRINTF("finished communicating setup_nper_rank_other_temp\n");
  }
  
  //common part of initialization
  void all_to_all_comm_t::common_setup_part1(temp_build_t& build)
  {
    //count the number of different ranks to send and receive from
    nranks_to=nranks_fr=0;
    for(int irank=0;irank<nranks;irank++)
      {
	if(build.nper_rank_to_temp[irank]!=0)
	  nranks_to++;
	
	if(build.nper_rank_fr_temp[irank]!=0)
	  nranks_fr++;
      }
    inited=true;
    list_ranks_to.resize(nranks_to);
    list_ranks_fr.resize(nranks_fr);
    
    //store the true list of ranks to communicate with and the number of elements to send to each rank
    nper_rank_to.resize(nranks_to);
    nper_rank_fr.resize(nranks_fr);
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
    out_buf_off_per_rank.resize(nranks_to);
    in_buf_off_per_rank .resize(nranks_fr);
    build.out_buf_cur_per_rank.resize(nranks_to);
    build.in_buf_cur_per_rank .resize(nranks_fr);
    
    if(nranks_to!=0)
      {
	if(nranks_to!=0)
	  {
	    out_buf_off_per_rank[0]=build.out_buf_cur_per_rank[0]=0;
	    for(int ilist_rank_to=1;ilist_rank_to<nranks_to;ilist_rank_to++)
	      out_buf_off_per_rank[ilist_rank_to]=build.out_buf_cur_per_rank[ilist_rank_to]=
		out_buf_off_per_rank[ilist_rank_to-1]+nper_rank_to[ilist_rank_to-1];
	  }
	if(nranks_fr!=0)
	  {
	    in_buf_off_per_rank[0]=build.in_buf_cur_per_rank[0]=0;
	    for(int ilist_rank_fr=1;ilist_rank_fr<nranks_fr;ilist_rank_fr++)
	      in_buf_off_per_rank[ilist_rank_fr]=build.in_buf_cur_per_rank[ilist_rank_fr]=
		in_buf_off_per_rank[ilist_rank_fr-1]+nper_rank_fr[ilist_rank_fr-1];
	  }
      }
  }
  
  //explain to each rank how to sort out or fill data from the buff
  //"expl" target the ranks to instruct
  //"note" means the ranks from which learning
  void all_to_all_comm_t::common_setup_part2(const int& nel_note,
					     int*& buf_note,
					     const int& nranks_note,
					     const std::vector<int>& list_ranks_note,
					     const std::vector<int>& buf_note_off_per_rank,
					     const std::vector<int>& nper_rank_note,
					     const std::vector<int>& buf_expl,
					     const int& nranks_expl,
					     const std::vector<int>& list_ranks_expl,
					     const std::vector<int>& buf_expl_off_per_rank,
					     const std::vector<int>& nper_rank_expl)
  {
    buf_note=nissa_malloc("buf_note",nel_note,int);
    
    MPI_Request req_list[nranks_note+nranks_expl];
    int ireq=0;
    for(int irank_expl=0;irank_expl<nranks_expl;irank_expl++)
      MPI_Isend(&buf_expl[buf_expl_off_per_rank[irank_expl]],
		nper_rank_expl[irank_expl],
		MPI_INT,
		list_ranks_expl[irank_expl],
		909,
		MPI_COMM_WORLD,
		&req_list[ireq++]);
    for(int irank_note=0;irank_note<nranks_note;irank_note++)
	  MPI_Irecv(&buf_note[buf_note_off_per_rank[irank_note]],
		    nper_rank_note[irank_note],
		    MPI_INT,
		    list_ranks_note[irank_note],
		    909,
		    MPI_COMM_WORLD,
		    &req_list[ireq++]);
    
    if(ireq!=nranks_note+nranks_expl)
      CRASH("expected %d request, obtained %d",nranks_note+nranks_expl,ireq);
    MPI_Waitall(ireq,req_list,MPI_STATUS_IGNORE);
    
    //check
    int max_nel_in=0;
    for(int iel_in=0;iel_in<nel_in;iel_in++)
      max_nel_in=std::max(max_nel_in,in_buf_dest[iel_in]);
    
    std::vector<int> in_buf_dest_check(max_nel_in+1,0);
    for(int iel_in=0;iel_in<nel_in;iel_in++)
      {
	int idest=in_buf_dest[iel_in];
	
	//if(idest<0 or idest>=nel_in) CRASH("in_buf_dest[%d] point to %d not in the range [0,nel_in=%d)",iel_in,idest,nel_in);
	if(in_buf_dest_check[idest]++==1)
	  CRASH("multiple assignement of %d",idest);
      }
  }
  
  //build knowing where to send
  all_to_all_comm_t::all_to_all_comm_t(const all_to_all_scattering_list_t& sl)
  {
    setup_knowing_where_to_send(sl);
  }
  
  void all_to_all_comm_t::setup_knowing_where_to_send(const all_to_all_scattering_list_t& sl)
  {
    //count the number of elements to send
    temp_build_t build;
    nel_out=sl.size();
    VERBOSITY_LV3_MASTER_PRINTF("nel to be scattered out: %d\n",nel_out);
    
    //count how many elements to send to each rank
    for(int& i : build.nper_rank_to_temp)
      i=0;
    
    for(int iel_out=0;iel_out<nel_out;iel_out++)
      {
	const int rank_to=
	  sl[iel_out].second%nranks;
	
	if(rank_to>=nranks or rank_to<0)
	  CRASH("destination rank %d does not exist!",rank_to);
	
	build.nper_rank_to_temp[rank_to]++;
      }
    
    //count how many elements to receive to each rank
    setup_nper_rank_other_temp(build.nper_rank_fr_temp,build.nper_rank_to_temp);
    common_setup_part1(build);
    out_buf_source=nissa_malloc("out_buf_source",nel_out,int);
    
    //save where to store and fill where to load each element in the temporary buffer
    std::vector<int> in_buf_dest_expl(nel_out);
    for(int iel_out=0;iel_out<nel_out;iel_out++)
      {
	const int rank_iel_to=sl[iel_out].second;
	
	const int iel_to=rank_iel_to/nranks;
	
	const int rank_to=rank_iel_to-nranks*iel_to;
	
	const int ilist_rank_to=build.rank_to_map_list_ranks_to[rank_to];
	
	const int ipos=build.out_buf_cur_per_rank[ilist_rank_to]++;
	
	out_buf_source[ipos]=sl[iel_out].first;
	
	in_buf_dest_expl[ipos]=iel_to;
      }
    
    common_setup_part2(nel_in,
		       in_buf_dest,
		       nranks_fr,
		       list_ranks_fr,
		       in_buf_off_per_rank,
		       nper_rank_fr,
		       in_buf_dest_expl,
		       nranks_to,
		       list_ranks_to,
		       out_buf_off_per_rank,
		       nper_rank_to);
  }
  
  //build knowing where to send
  all_to_all_comm_t::all_to_all_comm_t(const all_to_all_gathering_list_t& gl)
  {
    setup_knowing_what_to_ask(gl);
  }
  
  void all_to_all_comm_t::setup_knowing_what_to_ask(const all_to_all_gathering_list_t& gl)
  {
    temp_build_t build;
    nel_in=gl.size();
    
    //count how many elements to send to each rank
    for(int& i : build.nper_rank_fr_temp)
      i=0;
    
    for(auto it=gl.begin();it!=gl.end();it++)
      {
	int rank_fr=it->first%nranks;
	if(rank_fr>=nranks or rank_fr<0) CRASH("source rank %d does not exist!",rank_fr);
	build.nper_rank_fr_temp[rank_fr]++;
      }
    
    //count how many elements to receive from each rank
    setup_nper_rank_other_temp(build.nper_rank_to_temp,build.nper_rank_fr_temp);
    common_setup_part1(build);
    in_buf_dest=nissa_malloc("in_buf_dest",nel_in,int);
    
    //save the explenations to each rank on how to fill outbuffers
    std::vector<int> out_buf_source_expl(nel_in);
    for(auto it=gl.begin();it!=gl.end();it++)
      {
	const int rank_iel_fr=it->first;
	
	const int iel_fr=rank_iel_fr/nranks;
	
	const int rank_fr=rank_iel_fr-iel_fr*nranks;
	
	const int ilist_rank_fr=build.rank_fr_map_list_ranks_fr[rank_fr];
	
	const int ipos=build.in_buf_cur_per_rank[ilist_rank_fr]++;
	
	in_buf_dest[ipos]=it->second;
	  out_buf_source_expl[ipos]=iel_fr;
	}
    
    common_setup_part2(nel_out,
		       out_buf_source,
		       nranks_to,
		       list_ranks_to,
		       out_buf_off_per_rank,
		       nper_rank_to,
		       out_buf_source_expl,
		       nranks_fr,
		       list_ranks_fr,
		       in_buf_off_per_rank,
		       nper_rank_fr);
  }
  
  //perform the remapping
  void all_to_all_comm_t::communicate(void* out,
				      void* in,
				      const size_t& bps,
				      const int& tag) const
  {
    //allocate a buffer where to repack data
    char* out_buf=
      memoryManager<defaultMemorySpace>()->provide<char>(nel_out*bps);
    
    char* in_buf=
      memoryManager<defaultMemorySpace>()->provide<char>(nel_in*bps);
    
    const int* source=
      this->out_buf_source;
    
    //copy data on the out-going buffer
    PAR(0,nel_out,
	CAPTURE(bps,
		out_buf,
		source,
		in),
		iel_out,
		{
		  memcpy(out_buf+iel_out*bps,(char*)in+source[iel_out]*bps,bps);
		});
    
    /// Returns the argument after checking that does not exceed the max count
    const auto check_not_above_max_count=
      [](const size_t& n)
      {
	if(n>std::numeric_limits<int>::max())
	  CRASH("trying to send or receive %zu elements, max value is %d",n,std::numeric_limits<int>::max());
	
	return n;
      };
    
    MPI_Request req_list[nranks_to+nranks_fr];
    int ireq=0;
    int irank_fr_this=nranks_fr;
    for(int irank_fr=0;irank_fr<nranks_fr;irank_fr++)
      if(list_ranks_fr[irank_fr]!=rank)
	{
	  MPI_Irecv(in_buf+in_buf_off_per_rank[irank_fr]*bps,
		    check_not_above_max_count(nper_rank_fr[irank_fr]*bps),
		    MPI_CHAR,
		    list_ranks_fr[irank_fr],
		    909,
		    MPI_COMM_WORLD,&req_list[ireq++]);
	  // MASTER_PRINTF("Going to receive from rank %d, nreq: %d\n",list_ranks_fr[irank_fr],ireq);
	}
      else
	irank_fr_this=irank_fr;
    
    int irank_to_this=nranks_to;
    for(int irank_to=0;irank_to<nranks_to;irank_to++)
      if(list_ranks_to[irank_to]!=rank)
	{
	  MPI_Isend(out_buf+out_buf_off_per_rank[irank_to]*bps,
		    check_not_above_max_count(nper_rank_to[irank_to]*bps),
		    MPI_CHAR,
		    list_ranks_to[irank_to],
		    909,
		    MPI_COMM_WORLD,&req_list[ireq++]);
	  // MASTER_PRINTF("Going to send to rank %d, nreq: %d\n",list_ranks_to[irank_to],ireq);
	}
      else
	irank_to_this=irank_to;
    
    if(irank_fr_this!=nranks_fr and irank_to_this!=nranks_to)
      {
	if(nper_rank_to[irank_to_this]!=nper_rank_fr[irank_fr_this])
	  CRASH("unmatched what to copy locally");
	else
	  PAR(0,nper_rank_to[irank_fr_this],
	      CAPTURE(bps,
		      lhs=in_buf+in_buf_off_per_rank[irank_fr_this]*bps,
		      rhs=out_buf+out_buf_off_per_rank[irank_to_this]*bps),
	      i,
	      {
		memcpy(lhs+i*bps,rhs+i*bps,bps);
	      });
      }
    
    // MASTER_PRINTF("waiting for %d reqs\n",ireq);
    MPI_Waitall(ireq,req_list,MPI_STATUS_IGNORE);
    
    const int* dest=in_buf_dest;
    
    PAR(0,nel_out,
	CAPTURE(bps,
		in_buf,
		dest,
		out),
	iel_in,
	{
	  memcpy((char*)out+dest[iel_in]*bps,in_buf+iel_in*bps,bps);
	});
    
    memoryManager<defaultMemorySpace>()->release(out_buf);
    
    memoryManager<defaultMemorySpace>()->release(in_buf);
  }
    
    //add links to the buffer of the conf if needed
  int all_to_all_gathering_list_t::add_conf_link_for_paths(const Coords& g,
							   const int& mu)
  {
    //find rank and local position
    const auto[ivol,irank]=
      getLoclxAndRankOfCoords(g);
    
    const int ilink_asked=NDIM*ivol+mu;
    
    //if it is local, return local position
    if(irank==rank) return ilink_asked;
    else
      {
	const int irank_link_asked=ilink_asked*nranks+irank;
	
	//if it is non local search it in the list of to-be-gathered
	all_to_all_gathering_list_t::iterator it=this->find(irank_link_asked);
	
	//if it is already in the list, return its position
	if(it!=this->end()) return it->second;
	else
	  {
	    //otherwise add it to the list of to-be-gathered
	    int nel_gathered=this->size();
	    int igathered=NDIM*locVol+nel_gathered;
	    (*this)[irank_link_asked]=igathered;
	    
	    return igathered;
	  }
      }
  }
}
