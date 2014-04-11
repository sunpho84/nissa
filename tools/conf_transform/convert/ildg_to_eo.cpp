#include "nissa.hpp"

int T,L;

using namespace nissa;

void conf_convert(char *outpath,char *inpath)
{
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  read_ildg_gauge_conf(conf,inpath);
  
  FILE *fout=NULL;
  if(rank==0) fout=fopen(outpath,"w");  
  
  {
    coords temp;
    if(!little_endian) uint32s_to_uint32s_changing_endianness((uint32_t*)temp,(uint32_t*)glb_size,4);
    else                memcpy(temp,glb_size,sizeof(coords));
    if(rank==0)
      {
	int nw=fwrite(temp,sizeof(coords),1,fout);
	if(nw!=1) crash("did not success in writing");
      }
  }
  
  {
      double plaq=global_plaquette_lx_conf(conf)*3;
      if(!little_endian) doubles_to_doubles_changing_endianness(&plaq,&plaq,1);
      if(rank==0)
	{
	  int nw=fwrite(&plaq,sizeof(double),1,fout);
	  if(nw!=1) crash("did not success in writing");
	}
  }
   
  if(rank==0) fclose(fout);
  
  {
    if(!little_endian) doubles_to_doubles_changing_endianness((double*)conf,(double*)conf,loc_vol*4*18);
      char *buf=nissa_malloc("buf",loc_vol*sizeof(quad_su3),char);
      int ibuf=0;
      int remap[4]={0,3,2,1};
      coords g;
      int iscan=0;
      for(g[0]=0;g[0]<T;g[0]++)
	for(g[3]=0;g[3]<L;g[3]++)
	  for(g[2]=0;g[2]<L;g[2]++)
	    for(g[1]=0;g[1]<L;g[1]++)
	      if((g[0]+g[1]+g[2]+g[3])%2)
	        {
		  int dest_rank=iscan/loc_vol;
		  
		  int iw,iw_rank;
		  get_loclx_and_rank_of_coord(&iw,&iw_rank,g);
		  for(int mu=0;mu<4;mu++)
		    {
		      int nu=remap[mu];
		      int iz,iz_rank;
		      coords r;
		      memcpy(r,g,sizeof(coords));
		      r[nu]=(g[nu]+glb_size[nu]-1)%glb_size[nu];
		      get_loclx_and_rank_of_coord(&iz,&iz_rank,r);
		      
		      if(rank==iw_rank)
			{
			  if(dest_rank==rank)
			    {
			      memcpy(buf+ibuf,conf[iw][nu],sizeof(su3));
			      ibuf+=sizeof(su3);
			    }
			  else MPI_Send(conf[iw][remap[mu]],sizeof(su3),MPI_CHAR,dest_rank,0,MPI_COMM_WORLD);
			}
		      else
			if(rank==dest_rank)
			  {
			    MPI_Recv(buf+ibuf,sizeof(su3),MPI_CHAR,iw_rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			    ibuf+=sizeof(su3);
			  }

		      if(rank==iz_rank)
			{
			  if(dest_rank==rank)
			    {
			      memcpy(buf+ibuf,conf[iz][nu],sizeof(su3));
			      ibuf+=sizeof(su3);
			    }
			  else MPI_Send(conf[iz][remap[mu]],sizeof(su3),MPI_CHAR,dest_rank,0,MPI_COMM_WORLD);
			}
		      else
		      if(rank==dest_rank)
			{
			  MPI_Recv(buf+ibuf,sizeof(su3),MPI_CHAR,iz_rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			  ibuf+=sizeof(su3);
			}

		    }
		  iscan+=2;
		}
      if(ibuf!=(int)sizeof(quad_su3)*loc_vol) crash("did not arrive to the end: %d %d",ibuf,sizeof(quad_su3)*loc_vol);
  
      for(int irank=0;irank<nranks;irank++)
	{
	  if(rank==irank)
	    {
	      fout=fopen(outpath,"a");
	      
	      int nw=fwrite(buf,sizeof(quad_su3),loc_vol,fout);
	      /*
		NISSA_LOC_VOL_LOOP(ivol)
		for(int mu=0;mu<4;mu++)
		int nw=fprintf(fout,"%lg\n",((quad_su3*)buf)[ivol][mu][0][0][0]);
	      */
	      if(nw!=loc_vol) crash("did not success in writing");
	      
	      fflush(fout);
	      fclose(fout);
	    }
	  MPI_Barrier(MPI_COMM_WORLD);
	}
      nissa_free(buf);
  }

  nissa_free(conf);
}

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa(narg,arg);
  
  //if(nranks>1) crash("Cannot run in parallel");
  
  if(narg<2) crash("Use: %s input",arg[0]);
  
  open_input(arg[1]);
  
  // 1) Read information about the gauge conf
  
  //Read the volume
  read_str_int("L",&L);
  read_str_int("T",&T);
  
  //Init the MPI grid 
  init_grid(T,L);
  
  int N;
  read_str_int("NGaugeConf",&N);
  
  ///////////////////////////////////////////
  
  for(int i=0;i<N;i++)
  {
      char in[1024],out[1024];
      read_str(in,1024);
      read_str(out,1024);
      
      conf_convert(out,in);
  }
  
  close_input();
  
  ///////////////////////////////////////////
  
  close_nissa();
  
  return 0;
}
