#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <base/field.hpp>

namespace nissa
{
  std::vector<MPI_Request> startBufHaloNeighExchange(const int& divCoeff,
						     const size_t& bps)
  {
    /// Pending requests
    std::vector<MPI_Request> requests;
    
    for(int mu=0;mu<NDIM;mu++)
      if(is_dir_parallel[mu])
	for(int sendOri=0;sendOri<2;sendOri++)
	  {
	    const auto sendOrRecv=
	      [&mu,
	       &divCoeff,
	       &bps,
	       &requests,
	       messageTag=sendOri+2*mu]
	      (const char* oper,
	       const auto sendOrRecv,
	       auto* ptr,
	       const int& ori)
	      {
		const size_t offset=(bord_offset[mu]+bord_volh*ori)*bps/divCoeff;
		const int neighRank=rank_neigh[ori][mu];
		const size_t messageLength=bord_dir_vol[mu]*bps/divCoeff;
		// printf("rank %d %s ori %d dir %d, corresponding rank: %d, tag: %d length: %zu\n"
		//        ,rank,oper,ori,mu,neighRank,messageTag,messageLength);
		
		MPI_Request request;
		
		sendOrRecv(ptr+offset,messageLength,MPI_CHAR,neighRank,
			   messageTag,cart_comm,&request);
		
		requests.push_back(request);
	      };
	    
	    const int recvOri=1-sendOri;
	    
	    sendOrRecv("send",MPI_Isend,send_buf,sendOri);
	    sendOrRecv("recv",MPI_Irecv,recv_buf,recvOri);
	  }
    
    return requests;
  }
  
  /// Start the communications of buffer interpreted as edges
  std::vector<MPI_Request> startBufEdgesNeighExchange(const int& divCoeff,
						      const size_t& bps)
  {
    /// Pending requests
    std::vector<MPI_Request> requests;
    
    for(int iEdge=0;iEdge<nEdges;iEdge++)
      for(int sendOri1=0;sendOri1<2;sendOri1++)
	for(int sendOri2=0;sendOri2<2;sendOri2++)
	  {
	    const auto sendOrRecv=
	      [&iEdge,
	       &divCoeff,
	       &bps,
	       &requests,
	       messageTag=sendOri2+2*(sendOri1+2*iEdge)]
	       (const char* oper,
		const auto sendOrRecv,
		auto* ptr,
		const int& ori1,
		const int& ori2)
		  {
		    if(isEdgeParallel[iEdge])
		      {
			const size_t offset=(edge_offset[iEdge]+edge_vol*(ori2+2*ori1)/4)*bps/divCoeff;
			const int neighRank=rank_edge_neigh[ori1][ori2][iEdge];
			const size_t messageLength=edge_dir_vol[iEdge]*bps/divCoeff;
			
			// const auto [mu,nu]=edge_dirs[iEdge];
			// printf("rank %d %s ori %d,%d edge %d dir %d,%d, corresponding rank: %d, tag: %d length: %zu\n"
			//        ,rank,oper,ori1,ori2,iEdge,mu,nu,neighRank,messageTag,messageLength);
			
			MPI_Request request;
			
			sendOrRecv(ptr+offset,messageLength,MPI_CHAR,neighRank,
				   messageTag,cart_comm,&request);
			
			requests.push_back(request);
		      }
		  };
	    
	    const int recvOri1=1-sendOri1;
	    const int recvOri2=1-sendOri2;
	    
	    sendOrRecv("send",MPI_Isend,send_buf,sendOri1,sendOri2);
	    sendOrRecv("recv",MPI_Irecv,recv_buf,recvOri1,recvOri2);
	  }
    
    return requests;
  }

}
