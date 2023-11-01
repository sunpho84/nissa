#ifndef _ALL_TO_ALL_HPP
#define _ALL_TO_ALL_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <algorithm>
#include <map>
#include <vector>

#include <expr/comps.hpp>
#include <expr/dynamicCompsProvider.hpp>
#include <expr/dynamicTens.hpp>
#include <expr/mirroredNode.hpp>

namespace nissa
{
  template <typename Comps,
	    typename Fund,
	    bool IsRef=false>
  using MirroredTens=
    MirroredNode<Comps,DynamicTens<Comps,Fund,MemoryType::CPU,IsRef>>;
  
  /// All to all communicators
  template <typename CFr,
	    typename CTo,
	    bool IsRef=false>
  struct AllToAllComm
  {
    bool inited;
    
    struct BufComp :
      BaseComp<BufComp,int64_t,0>
    {
      using Base=BaseComp<BufComp,int64_t,0>;
      using Base ::Base;
    };
    
    struct Rank :
      BaseComp<Rank,int64_t,0>
    {
      using Base=BaseComp<Rank,int64_t,0>;
      using Base ::Base;
    };
    
    AllToAllComm() :
      inited{false}
    {
    }
    
    MirroredTens<CompsList<CFr>,BufComp,IsRef> inBufDest;
    
    MirroredTens<CompsList<BufComp>,CTo,IsRef> outBufDest;
    
    template <typename F>
    void init(const CFr& nCfr,
	      F&& f)
    {
      if(inited)
	crash("Cannot init twice");
      
      std::vector<CompsList<Rank,CTo>> sl;
      for(CFr cFr=0;cFr<nCfr;cFr++)
	{
	  const CompsList<Rank,CTo> rcTo=f(cFr);
	  
	  const Rank& rTo=std::get<Rank>(rcTo);
	  if(rTo>=nranks or rTo<0)
	    crash("destination rank %d does not exist!",rTo());
	  
	  sl.push_back(rcTo);
	}
      
      DynamicTens<CompsList<Rank>,int64_t,MemoryType::CPU> nperRankToTemp((Rank(nranks)));
      
      DynamicTens<CompsList<Rank>,int64_t,MemoryType::CPU> nperRankFromTemp((Rank(nranks)));

      int *out_buf_cur_per_rank,*in_buf_cur_per_rank;
      std::map<int,int> rank_to_map_list_ranks_to,rank_fr_map_list_ranks_fr;

      nel_out=sl.size();
    verbosity_lv3_master_printf("nel to be scattered out: %d\n",nel_out);
    
    //count how many elements to send to each rank
    vector_reset(build.nper_rank_to_temp);
    for(int iel_out=0;iel_out<nel_out;iel_out++)
      {
	int rank_to=sl[iel_out].second%nranks;
	if(rank_to>=nranks or rank_to<0) crash("destination rank %d does not exist!",rank_to);
	build.nper_rank_to_temp[rank_to]++;
      }
    
    //count how many elements to receive to each rank
    //send and receive the number of elements to communicate with
    for(int delta_rank=0;delta_rank<nranks;delta_rank++)
      {
	int dest_rank=(rank+nranks+delta_rank)%nranks;
	int recv_rank=(rank+nranks-delta_rank)%nranks;
	MPI_Sendrecv(nper_rank_temp+dest_rank,1,MPI_INT,dest_rank,0,
		     nper_rank_other_temp+recv_rank,1,MPI_INT,recv_rank,0,
		     MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
    
    verbosity_lv3_master_printf("finished communicating setup_nper_rank_other_temp\n");

    //count the number of different ranks to send and receive from
    int nranks_to_loc=0,nranks_fr_loc=0; //local: structure could be global!
    for(int irank=0;irank<nranks;irank++)
      {
	if(build.nper_rank_to_temp[irank]!=0) nranks_to_loc++;
	if(build.nper_rank_fr_temp[irank]!=0) nranks_fr_loc++;
      }
    inited=true;
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
    out_buf_source=nissa_malloc("out_buf_source",nel_out,int);
    
    //save where to store and fill where to load each element in the temporary buffer
    int *in_buf_dest_expl=nissa_malloc("in_buf_dest_expl",nel_out,int);
    
    for(int iel_out=0;iel_out<nel_out;iel_out++)
      {
	int rank_iel_to=sl[iel_out].second;
	int iel_to=rank_iel_to/nranks,rank_to=rank_iel_to-nranks*iel_to;
	int ilist_rank_to=build.rank_to_map_list_ranks_to[rank_to];
	int ipos=build.out_buf_cur_per_rank[ilist_rank_to]++;
	out_buf_source[ipos]=sl[iel_out].first;
	in_buf_dest_expl[ipos]=iel_to;
      }
    
    buf_note=nissa_malloc("buf_note",nel_note,int);
    
    MPI_Request req_list[nranks_note+nranks_expl];
    int ireq=0;
    for(int irank_expl=0;irank_expl<nranks_expl;irank_expl++)
      MPI_Isend(buf_expl+buf_expl_off_per_rank[irank_expl],nper_rank_expl[irank_expl],MPI_INT,
		list_ranks_expl[irank_expl],909,MPI_COMM_WORLD,&req_list[ireq++]);
    for(int irank_note=0;irank_note<nranks_note;irank_note++)
	  MPI_Irecv(buf_note+buf_note_off_per_rank[irank_note],nper_rank_note[irank_note],MPI_INT,
		    list_ranks_note[irank_note],909,MPI_COMM_WORLD,&req_list[ireq++]);
    if(ireq!=nranks_note+nranks_expl) crash("expected %d request, obtained %d",nranks_note+nranks_expl,ireq);
    MPI_Waitall(ireq,req_list,MPI_STATUS_IGNORE);
    
    //check
    int max_nel_in=0;
    for(int iel_in=0;iel_in<nel_in;iel_in++) max_nel_in=std::max(max_nel_in,in_buf_dest[iel_in]);
    
    int *in_buf_dest_check=nissa_malloc("in_buf_dest_check",max_nel_in+1,int);
    vector_reset(in_buf_dest_check);
    for(int iel_in=0;iel_in<nel_in;iel_in++)
      {
	int idest=in_buf_dest[iel_in];
	
	//if(idest<0 or idest>=nel_in) crash("in_buf_dest[%d] point to %d not in the range [0,nel_in=%d)",iel_in,idest,nel_in);
	if(in_buf_dest_check[idest]++==1) crash("multiple assignement of %d",idest);
      }
    
    nissa_free(in_buf_dest_check);
    nissa_free(in_buf_dest_expl);
      // setup_knowing_where_to_send(sl);
      
      inited=true;
    }
    
    // int nel_out{0},nel_in{0};
    // int nranks_fr{0},*list_ranks_fr{nullptr},*in_buf_dest{nullptr},*nper_rank_fr{nullptr},*in_buf_off_per_rank{nullptr};
    // int nranks_to{0},*list_ranks_to{nullptr},*out_buf_source{nullptr},*nper_rank_to{nullptr},*out_buf_off_per_rank{nullptr};
    
    // all_to_all_comm_t(const all_to_all_gathering_list_t &gl);
    // all_to_all_comm_t(const all_to_all_scattering_list_t &sl);
    // all_to_all_comm_t(all_to_all_comm_t &&oth)
    // {
    //   std::swap(inited,oth.inited);
    //   std::swap(nel_out,oth.nel_out);
    //   std::swap(nel_in,oth.nel_in);
    //   std::swap(nranks_fr,oth.nranks_fr);
    //   std::swap(list_ranks_fr,oth.list_ranks_fr);
    //   std::swap(in_buf_dest,oth.in_buf_dest);
    //   std::swap(nper_rank_fr,oth.nper_rank_fr);
    //   std::swap(in_buf_off_per_rank,oth.in_buf_off_per_rank);
    //   std::swap(nranks_to,oth.nranks_to);
    //   std::swap(list_ranks_to,oth.list_ranks_to);
    //   std::swap(out_buf_source,oth.out_buf_source);
    //   std::swap(nper_rank_to,oth.nper_rank_to);
    //   std::swap(out_buf_off_per_rank,oth.out_buf_off_per_rank);
    // }
    
    // ~all_to_all_comm_t(){destroy();}
    // void destroy()
    // {
    //   if(inited)
    // 	{
    // 	  inited=false;
	  
    // 	  nissa_free(list_ranks_to);
    // 	  nissa_free(list_ranks_fr);
    // 	  nissa_free(in_buf_dest);
    // 	  nissa_free(out_buf_source);
    // 	  nissa_free(nper_rank_fr);
    // 	  nissa_free(nper_rank_to);
    // 	  nissa_free(out_buf_off_per_rank);
    // 	  nissa_free(in_buf_off_per_rank);
    // 	}
    // }
    
    // void communicate(void *out,
    // 		     void *in,
    // 		     size_t bps,
    // 		     void *buf_out=NULL,
    // 		     void *buf_in=NULL,
    // 		     int tag=-1) const;
    
    // void setup_knowing_where_to_send(const all_to_all_scattering_list_t &sl);
    // void setup_knowing_what_to_ask(const all_to_all_gathering_list_t &gl);
    // void setup_nper_rank_other_temp(int *nper_rank_other_temp,int *nper_rank_temp);
    // void common_setup_part1(temp_build_t &build);
    // void common_setup_part2(int nel_note,int *&buf_note,int nranks_note,int *list_ranks_note,int *buf_note_off_per_rank,int *nper_rank_note,int *buf_expl,int nranks_expl,int *list_ranks_expl,int *buf_expl_off_per_rank,int *nper_rank_expl);
  };
}

#endif
