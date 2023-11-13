#ifndef _ALL_TO_ALL_HPP
#define _ALL_TO_ALL_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

/// \file expr/allToAll.hpp

#include <algorithm>
#include <map>
#include <vector>

#include <expr/comps.hpp>
#include <expr/dynamicCompsProvider.hpp>
#include <expr/dynamicTens.hpp>
#include <expr/mirroredNode.hpp>
#include <threads/threads.hpp>
#include <tuples/tupleCat.hpp>

namespace nissa
{
  /// Tensor to be used to store buffer source and dest
  template <typename Comps,
	    typename Fund,
	    bool IsRef=false>
  using MirroredTens=
    MirroredNode<Comps,DynamicTens<Comps,Fund,MemoryType::CPU,IsRef>>;
  
  /// All to all communicators
  template <DerivedFromComp CDst,
	    DerivedFromComp CSrc,
	    bool IsRef=false>
  struct AllToAllComm
  {
    /// Status: is the communicator inited?
    bool inited;
    
    /// Component to be used for the buffer
    struct BufComp :
      BaseComp<BufComp,int64_t,0>
    {
      using Base=BaseComp<BufComp,int64_t,0>;
      using Base::Base;
    };
    
    /// Default constructor
    AllToAllComm() :
      inited{false}
    {
    }
    
    /// Destination of the source, in the buffer
    MirroredTens<CompsList<CSrc>,BufComp,IsRef> outBufOfSrc;
    
    /// Destination of the buffer, in the final storage
    MirroredTens<CompsList<BufComp>,CDst,IsRef> dstOfInBuf;
    
    /// Keeps track of the number of elements to be sent to each given rank
    std::vector<std::pair<MpiRank,int64_t>> nSendToRank;
    
    /// Keeps track of the number of elements to be received from each given rank
    std::vector<std::pair<MpiRank,int64_t>> nRecvFrRank;
    
    /// Gets the destination size, which is equal to the in buffer size
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    CDst getNDst() const
    {
      return getInBufSize()();
    }
    
    /// Gets the incoming buffer size
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    constexpr BufComp getInBufSize() const
    {
      return dstOfInBuf.template getCompSize<BufComp>();
    }
    
    /// Gets the source size, which is equal to the out buffer size
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    CSrc getNSrc() const
    {
      return getOutBufSize()();
    }
    
    /// Gets the outgoing buffer size
    CUDA_HOST_AND_DEVICE INLINE_FUNCTION
    constexpr BufComp getOutBufSize() const
    {
      return outBufOfSrc.template getCompSize<CSrc>()();
    }
    
    /// Initializes the AllToAll communicator
    template <typename F>
    void init(const CSrc& nSrc,
	      F&& f)
    {
      if(inited)
	crash("Cannot init twice");
      inited=true;
      
      outBufOfSrc.allocate(nSrc);
      
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
      
      // for(MpiRank iRank=0;iRank<nranks;iRank++)
      // 	printf("AllToAll Rank %d transmitting to rank %d: %zu el\n",rank,iRank(),locSrcsGroupedByDstRank[iRank()].size());
      
      /// Progressive storage of the dst
      std::vector<CDst> buildDstOfInBuf;
      
      BufComp nOutBuf=0;
      for(int dRank=0;dRank<nranks;dRank++)
	{
	  /// Rank towards which to send
	  const int sendRank=
	    (rank+nranks+dRank)%nranks;
	  
	  for(const CSrc& locSrc : locSrcsGroupedByDstRank[sendRank])
	    {
	      // printf("AllToAll Rank %d filling %zu with nOutBuf %zu to be sent to rank %d\n",rank,locSrc(),nOutBuf(),sendRank);
	      outBufOfSrc[locSrc]=nOutBuf++;
	    }
	  
	  /// Rank from which to receive
	  const int recvRank=
	    (rank+nranks-dRank)%nranks;
	  
	  /// List of position in the ouput where to store the portion of the buffer relative to recvRank
	  const std::vector<CDst> dstOfBufFrRank=
	    mpiSendrecv(sendRank,remDstsGroupedByDstRank[sendRank],recvRank);
	  
	  // printf("AllToAll Rank %d dRank %d sent to rank %d value
	  // nRemDstOfRank=%zu received from rank %d value
	  // %zu\n",rank,dRank,sendRank,remDstsGroupedByDstRank[sendRank].size(),recvRank,dstOfBufFrRank.size());
	  
	  // Increment the list of object to pull out from the buffer
	  for(const CDst& bufDst : dstOfBufFrRank)
	    buildDstOfInBuf.emplace_back(bufDst);
	  
	  // If the node has to receive, mark the node in the list with the size attached
	  if(const size_t s=dstOfBufFrRank.size();s)
	    nRecvFrRank.emplace_back(recvRank,s);
	  
	  // If the node has to send, mark the node to the list with the size attached
	  if(const size_t s=remDstsGroupedByDstRank[sendRank].size();s)
	    nSendToRank.emplace_back(sendRank,s);
	}
      
      dstOfInBuf.allocate((BufComp)(int64_t)buildDstOfInBuf.size());
      for(BufComp bc=0;const CDst& cDst : buildDstOfInBuf)
	dstOfInBuf[bc++]=cDst;
      
      dstOfInBuf.updateDeviceCopy();
      outBufOfSrc.updateDeviceCopy();
      
      std::vector<BufComp> inBufOfDest(getNDst()(),-1);
      for(BufComp i=0;i<getInBufSize();i++)
	{
	  auto& p=inBufOfDest[dstOfInBuf[i]()];
	  if(p!=-1)
	    crash("On rank %d destination %ld is filled by %ld and at least %ld at the same time",
		  rank,(int64_t)dstOfInBuf[i](),(int64_t)p(),(int64_t)i());
	  p=i;
	}
    }
    
    /// Communicate the expression
    template <DerivedFromNode _DstExpr,
	      DerivedFromNode SrcExpr>
    void communicate(_DstExpr&& out,
		     const SrcExpr& in) const
    {
      using DstExpr=std::decay_t<_DstExpr>;
      
      static_assert(tupleHasType<typename DstExpr::Comps,CDst>,"Destination does not have the destination component");
      static_assert(tupleHasType<typename SrcExpr::Comps,CSrc>,"Source does not have the source component");
      
      constexpr MemoryType dstExecSpace=DstExpr::execSpace;
      constexpr MemoryType srcExecSpace=SrcExpr::execSpace;
      constexpr MemoryType execSpace=srcExecSpace;
      
      static_assert(dstExecSpace==srcExecSpace,"Needs to have same exec space for src and dest");
      
      using DstRedComps=TupleFilterAllTypes<typename DstExpr::Comps,std::tuple<CDst>>;
      using SrcRedComps=TupleFilterAllTypes<typename SrcExpr::Comps,std::tuple<CSrc>>;
      static_assert(tupleHaveTypes<DstRedComps,SrcRedComps>,"Components in the source must be the same apart from CDst and CSrc");
      
      using Fund=DstExpr::Fund;
      
      using BufComps=
	TupleCat<DstRedComps,CompsList<BufComp>>;
      
      DynamicTens<BufComps,Fund,execSpace> outBuf(std::tuple_cat(std::make_tuple(getOutBufSize()),in.getDynamicSizes()));
      
      const CSrc nSrc=in.template getCompSize<CSrc>();
      if(nSrc()!=getOutBufSize()())
	crash("Size of the in epxression %ld different from expected %ld\n",(int64_t)nSrc(),(int64_t)getOutBufSize()());
      
      PAR_ON_EXEC_SPACE(execSpace,
			0,
			nSrc,
			CAPTURE(TO_READ(in),
				outBufOfSrc=outBufOfSrc.template getRefForExecSpace<execSpace>(),
				TO_WRITE(outBuf)),
			iIn,
			{
			  outBuf(outBufOfSrc(iIn))=in(iIn);
			});
      
      decltype(auto) hostOutBuf=
	outBuf.template copyToMemorySpaceIfNeeded<MemoryType::CPU>();
      const auto mergedHostOutBuf=
	hostOutBuf.template mergeComps<DstRedComps>();
      
      DynamicTens<BufComps,Fund,MemoryType::CPU> hostInBuf(std::tuple_cat(std::make_tuple(getInBufSize()),out.getDynamicSizes()));
      auto mergedHostInBuf=
	hostInBuf.template mergeComps<DstRedComps>();
      
      const int nReqs=nSendToRank.size()+nRecvFrRank.size();
      MPI_Request reqs[nReqs];
      MPI_Request* req=reqs;
      const auto nDof=mergedHostOutBuf.template getCompSize<MergedComp<DstRedComps>>()();
      for(BufComp sendOffset=0;const auto& [dstRank,nEl] : nSendToRank)
	{
	  MPI_Isend(&mergedHostOutBuf(sendOffset,MergedComp<DstRedComps>(0)),
		    nEl*nDof*sizeof(Fund),
		    MPI_CHAR,dstRank(),0,MPI_COMM_WORLD,req++);
	  sendOffset+=nEl;
	}
      
      for(BufComp recvOffset=0;const auto& [rcvRank,nEl] : nRecvFrRank)
	{
	  MPI_Irecv(&mergedHostInBuf(recvOffset,MergedComp<DstRedComps>(0)),
		    nEl*nDof*sizeof(Fund),
		    MPI_CHAR,rcvRank(),0,MPI_COMM_WORLD,req++);
	  recvOffset+=nEl;
	}
      MPI_Waitall(nReqs,reqs,MPI_STATUS_IGNORE);
      
      decltype(auto) inBuf=
	hostInBuf.template copyToMemorySpaceIfNeeded<execSpace>();
      
      PAR_ON_EXEC_SPACE(execSpace,
			0,
			getInBufSize(),
			CAPTURE(TO_READ(inBuf),
				dstOfInBuf=dstOfInBuf.template getRefForExecSpace<execSpace>(),
				TO_WRITE(out)),
			iInBuf,
			{
			  const CDst iOut=dstOfInBuf(iInBuf);
			  
			  out(iOut)=inBuf(iInBuf);
			});
      
    }
  };
}

#endif
