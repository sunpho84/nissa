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
  void vector_gather(char *glb,char *loc,size_t bps,int dest_rank)
  {
    crash("reimplement");
    
    // if(dest_rank==rank)
    //   {
    // 	//copy local data
    // 	memcpy(glb+rank*locVol*bps,loc,locVol*bps);
	
    // 	//open incoming communications for non-local data
    // 	MPI_Request req[nranks-1];
    // 	int ireq=0;
    // 	for(int irank=0;irank<nranks;irank++)
    // 	  if(irank!=rank)
    // 	    MPI_Irecv(glb+irank*locVol*bps,locVol*bps,MPI_CHAR,irank,239+irank,MPI_COMM_WORLD,&(req[ireq++]));
    // 	//wait all incoming data
    // 	MPI_Waitall(ireq,req,MPI_STATUS_IGNORE);
	
    // 	//reorder data
    // 	int *ord=nissa_malloc("ord",glbVol,int);
    // 	Coords r;
    // 	for(r[0]=0;r[0]<nrank_dir[0];r[0]++)
    // 	  for(r[1]=0;r[1]<nrank_dir[1];r[1]++)
    // 	    for(r[2]=0;r[2]<nrank_dir[2];r[2]++)
    // 	      for(r[3]=0;r[3]<nrank_dir[3];r[3]++)
    // 		{
    // 		  int irank=rank_of_coord(r);
    // 		  Coords l;
    // 		  for(l[0]=0;l[0]<locSize[0];l[0]++)
    // 		    for(l[1]=0;l[1]<locSize[1];l[1]++)
    // 		      for(l[2]=0;l[2]<locSize[2];l[2]++)
    // 			for(l[3]=0;l[3]<locSize[3];l[3]++)
    // 			  {
    // 			    Coords g;
    // 			    for(int mu=0;mu<4;mu++) g[mu]=r[mu]*locSize[mu]+l[mu];
			    
    // 			    int ivol=locVol*irank+loclx_of_coord(l);
    // 			    int glb_site=glblx_of_coord(g);
			    
    // 			    ord[ivol]=glb_site;
    // 			  }
    // 		}
	
    // 	reorder_vector(glb,ord,glbVol,bps);
	
    // 	nissa_free(ord);
    //   }
    // else
    //   {
    // 	//send non-local data
    // 	MPI_Request req;
    // 	MPI_Isend(loc,locVol*bps,MPI_CHAR,dest_rank,239+rank,MPI_COMM_WORLD,&req);
    // 	MPI_Waitall(1,&req,MPI_STATUS_IGNORE);
    //   }
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
    if(glbSize[0]%2 || glbSize[1]%2) crash("Error, impossible to mirrorize if sites are odds");
    
    int TH=glbSize[0]/2;
    int LH=glbSize[1]/2;
    
    int x[4];
    for(x[0]=0;x[0]<=TH;x[0]++)
      for(x[1]=0;x[1]<=LH;x[1]++)
	for(x[2]=0;x[2]<=LH;x[2]++)
	  for(x[3]=0;x[3]<=LH;x[3]++)
	    {
	      //find ipercubic mirrored partners
	      int ivol[8];
	      for(int imirr=0;imirr<8;imirr++)
		{
		  Coords xmirr;
		  for(int mu=0;mu<4;mu++)
		    xmirr[mu]=(imirr & (1<<mu)) ? (glbSize[mu]-x[mu])%glbSize[mu] : x[mu];
		  ivol[imirr]=glblxOfCoord(xmirr);
		}
	      
	      //average
	      average_list_of_gathered_vector_sites(vec,ivol,8,dps);
	    }
  }
  
  //symmetrize an already gathered vector
  void gathered_vector_cubic_symmetrize(double *vec,int dps)
  {
    if(glbSize[0]%2 || glbSize[1]%2) crash("Error, impossible to symmetrize if sites are odds");
    
    int TH=glbSize[0]/2;
    int LH=glbSize[1]/2;
    
    int perm[6][3]={{1,2,3},{2,3,1},{3,1,2},{1,3,2},{3,2,1},{2,1,3}};
    
    int x[4];
    for(x[0]=0;x[0]<=TH;x[0]++)
      for(x[1]=0;x[1]<=LH;x[1]++)
	for(x[2]=0;x[2]<=x[1];x[2]++)
	  for(x[3]=0;x[3]<=x[2];x[3]++)
	    {
	      //find cubic partners
	      int ivol[6];
	      for(int iperm=0;iperm<6;iperm++)
		ivol[iperm]=glblxOfCoordList(x[0],x[perm[iperm][0]],x[perm[iperm][1]],x[perm[iperm][2]]);
	      
	      //average
	      average_list_of_gathered_vector_sites(vec,ivol,6,dps);
	    }
  }
}
