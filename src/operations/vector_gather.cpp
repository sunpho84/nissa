#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"

namespace nissa
{
  //gather the whole field on a single rank, reordering data
  void vector_gather(char *glb,char *loc,const int64_t& bps,int dest_rank)
  {
    if(dest_rank==rank)
      {
	//copy local data
	memcpy(glb+(rank*locVol*bps).nastyConvert(),loc,(locVol*bps).nastyConvert());
	
	//open incoming communications for non-local data
	MPI_Request req[nranks-1];
	int ireq=0;
	for(int irank=0;irank<nranks;irank++)
	  if(irank!=rank)
	    MPI_Irecv(glb+(irank*locVol*bps).nastyConvert(),(locVol*bps).nastyConvert(),MPI_CHAR,irank,239+irank,MPI_COMM_WORLD,&(req[ireq++]));
	//wait all incoming data
	MPI_Waitall(ireq,req,MPI_STATUS_IGNORE);
	
	//reorder data
	int *ord=nissa_malloc("ord",glbVol.nastyConvert(),int);
	RankCoords r;
	for(Rank irank=0;irank<nranks;irank++)
	  for(LocLxSite ivol=0;ivol<locVol;ivol++)
	    {
	      GlbCoords g;
	      FOR_ALL_DIRECTIONS(mu)
		g(mu)=r(mu)()*locSize(mu)()+locCoordOfLoclx(ivol,mu)();
			    
	      const GlbLxSite& glb_site=glblx_of_coord(g);
			    
	      ord[ivol.nastyConvert()]=glb_site.nastyConvert();
	    }
	
	reorder_vector(glb,ord,glbVol.nastyConvert(),bps);
	
	nissa_free(ord);
      }
    else
      {
	//send non-local data
	MPI_Request req;
	MPI_Isend(loc,(locVol*bps).nastyConvert(),MPI_CHAR,dest_rank,239+rank,MPI_COMM_WORLD,&req);
	MPI_Waitall(1,&req,MPI_STATUS_IGNORE);
      }
  }
  
  //average over all the passed sites
  void average_list_of_gathered_vector_sites(double *vec,int *sites,int nsites,int dps)
  {
    double buf[dps];
    memcpy(buf,vec+dps*sites[0],sizeof(double)*dps);
    for(int id=0;id<dps;id++)
      {	
	for(int isite=1;isite<nsites;isite++)
	  buf[id]+=vec[dps*sites[isite]+id];
	buf[id]/=nsites;
      }
    
    //copy the symmetrized buffer over all the different partners
    for(int isite=0;isite<nsites;isite++)
      memcpy(vec+dps*sites[isite],buf,sizeof(double)*dps);
  }
  
  //ipercubicly mirrorize an already gathered vector
  void gathered_vector_mirrorize(double *vec,int dps)
  {
    crash("reimplement");
    // if(glbSize[0]%2 or glbSize[1]%2) crash("Error, impossible to mirrorize if sites are odds");
    
    // int TH=glbSize[0]/2;
    // int LH=glbSize[1]/2;
    
    // int x[4];
    // for(x[0]=0;x[0]<=TH;x[0]++)
    //   for(x[1]=0;x[1]<=LH;x[1]++)
    // 	for(x[2]=0;x[2]<=LH;x[2]++)
    // 	  for(x[3]=0;x[3]<=LH;x[3]++)
    // 	    {
    // 	      //find ipercubic mirrored partners
    // 	      int ivol[8];
    // 	      for(int imirr=0;imirr<8;imirr++)
    // 		{
    // 		  int xmirr[4];
    // 		  for(int mu=0;mu<4;mu++)
    // 		    xmirr[mu]=(imirr & (1<<mu)) ? (glbSize[mu]-x[mu])%glbSize[mu] : x[mu];
    // 		  ivol[imirr]=glblx_of_coord(xmirr);
    // 		}
	      
    // 	      //average
    // 	      average_list_of_gathered_vector_sites(vec,ivol,8,dps);
    // 	    }
  }
  
  //symmetrize an already gathered vector
  void gathered_vector_cubic_symmetrize(double *vec,int dps)
  {
    crash("reimplement");
    // if(glbSize[0]%2 || glbSize[1]%2) crash("Error, impossible to symmetrize if sites are odds");
    
    // int TH=glbSize[0]/2;
    // int LH=glbSize[1]/2;
    
    // int perm[6][3]={{1,2,3},{2,3,1},{3,1,2},{1,3,2},{3,2,1},{2,1,3}};
    
    // int x[4];
    // for(x[0]=0;x[0]<=TH;x[0]++)
    //   for(x[1]=0;x[1]<=LH;x[1]++)
    // 	for(x[2]=0;x[2]<=x[1];x[2]++)
    // 	  for(x[3]=0;x[3]<=x[2];x[3]++)
    // 	    {
    // 	      //find cubic partners
    // 	      int ivol[6];
    // 	      for(int iperm=0;iperm<6;iperm++)
    // 		ivol[iperm]=glblx_of_coord_list(x[0],x[perm[iperm][0]],x[perm[iperm][1]],x[perm[iperm][2]]);
	      
    // 	      //average
    // 	      average_list_of_gathered_vector_sites(vec,ivol,6,dps);
    // 	    }
  }
}
