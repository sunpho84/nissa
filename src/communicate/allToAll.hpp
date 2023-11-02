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
  template <typename CSrc,
	    typename CDst,
	    bool IsRef=false>
  struct AllToAllComm
  {
    bool inited;
    
    struct BufComp :
      BaseComp<BufComp,int64_t,0>
    {
      using Base=BaseComp<BufComp,int64_t,0>;
      using Base::Base;
    };
    
    AllToAllComm() :
      inited{false}
    {
    }
    
    MirroredTens<CompsList<CSrc>,BufComp,IsRef> outBufOfSrc;
    
    MirroredTens<CompsList<BufComp>,CDst,IsRef> dstOfInBuf;
    
    std::vector<std::pair<MpiRank,size_t>> nSendToRank;
    
    std::vector<std::pair<MpiRank,size_t>> nRecvFrRank;
    
    template <typename F>
    void init(const CSrc& nSrc,
	      const CDst& nDst,
	      F&& f)
    {
      if(inited)
	crash("Cannot init twice");
      inited=true;
      
      outBufOfSrc.allocate(nSrc);
      dstOfInBuf.allocate(BufComp(nDst()));
      
      /// For each destination rank, list all local sources
      std::vector<std::vector<CSrc>> locSrcsGroupedByDstRank(nranks);
      
      /// For each destination rank, list all remote destinations
      std::vector<std::vector<CDst>> remDstsGroupedByDstRank(nranks);
      for(CSrc locSrc=0;locSrc<nSrc;locSrc++)
	{
	  const auto [dstRank,remDst]=f(locSrc);
	  
	  if(dstRank>=nranks or dstRank<0)
	    crash("destination rank %d does not exist!",dstRank);
	  
	  locSrcsGroupedByDstRank[dstRank()].emplace_back(locSrc);
	  remDstsGroupedByDstRank[dstRank()].emplace_back(remDst);
	}
      
      for(MpiRank iRank=0;iRank<nranks;iRank++)
	printf("AllToAll Rank %d transmitting to rank %d: %zu el\n",rank,iRank(),locSrcsGroupedByDstRank[iRank()].size());
      
      BufComp nInBuf=0,nOutBuf=0;
      for(int dRank=0;dRank<nranks;dRank++)
	{
	  const int sendRank=(rank+nranks+dRank)%nranks;
	  const int recvRank=(rank+nranks-dRank)%nranks;
	  
	  const size_t nRemDstOfRank=
	    remDstsGroupedByDstRank[sendRank].size();
	  
	  size_t nRcvFrRank;
	  MPI_Sendrecv(&nRemDstOfRank,1,MPI_UINT64_T,sendRank,0,
		       &nRcvFrRank,   1,MPI_UINT64_T,recvRank,0,
		       MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  
	  printf("AllToAll Rank %d dRank %d sent to rank %d value nRemDstOfRank=%zu received from rank %d value %zu\n",rank,dRank,sendRank,nRemDstOfRank,recvRank,nRcvFrRank);
	  
	  std::vector<CDst> dstOfBufFrRank(nRcvFrRank);
	  MPI_Sendrecv(&remDstsGroupedByDstRank[0],nRemDstOfRank*sizeof(CDst),MPI_CHAR,sendRank,0,
		       &dstOfBufFrRank[0],nRcvFrRank*sizeof(CDst),MPI_CHAR,recvRank,0,
		       MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  
	  for(const CDst& bufDst : dstOfBufFrRank)
	    dstOfInBuf[nInBuf++]=bufDst;
	  
	  for(const CSrc& locSrc : locSrcsGroupedByDstRank[sendRank])
	    outBufOfSrc[locSrc]=nOutBuf++;
	  
	  if(const size_t s=dstOfBufFrRank.size();s)
	    nRecvFrRank.emplace_back(MpiRank(recvRank),s);
	  
	  if(const size_t s=remDstsGroupedByDstRank.size();s)
	    nSendToRank.emplace_back(MpiRank(sendRank),s);
	}
      
      dstOfInBuf.updateDeviceCopy();
      outBufOfSrc.updateDeviceCopy();
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
