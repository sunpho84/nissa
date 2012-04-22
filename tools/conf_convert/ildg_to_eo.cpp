#include "nissa.h"

int T,L;

void conf_convert(char *outpath,char *inpath)
{
  quad_su3 *conf=nissa_malloc("conf",loc_vol,quad_su3);
  read_ildg_gauge_conf(conf,inpath);
  
  FILE *fout=fopen(outpath,"wb");
  
  {
      coords temp;
      if(!little_endian) ints_to_ints_changing_endianess(temp,glb_size,4);
      else                memcpy(temp,glb_size,sizeof(coords));
      int nw=fwrite(temp,sizeof(coords),1,fout);
      if(nw!=1) crash("did not success in writing");
  }
  
  {
      double plaq=global_plaquette_lx_conf(conf)*3;
      if(!little_endian) doubles_to_doubles_changing_endianess(&plaq,&plaq,1);
      int nw=fwrite(&plaq,sizeof(double),1,fout);
      if(nw!=1) crash("did not success in writing");
  }
  
  {
      if(!little_endian) doubles_to_doubles_changing_endianess((double*)conf,(double*)conf,glb_vol*4*18);
      char *buf=nissa_malloc("buf",glb_vol*sizeof(quad_su3),char);
      int ibuf=0;
      int remap[4]={0,3,2,1};
      for(size_t t=0;t<T;t++)
	for(size_t z=0;z<L;z++)
	  for(size_t y=0;y<L;y++)
	    for(size_t x=0;x<L;x++)
	      if((x+y+z+t)%2)
	        {
		  int iw=loclx_of_coord_list(t,x,y,z);
		  for(int mu=0;mu<4;mu++)
		    {
		      int iz=loclx_neighdw[iw][remap[mu]];
		      
		      memcpy(buf+ibuf,conf[iw][remap[mu]],sizeof(su3));
		      ibuf+=sizeof(su3);
		      memcpy(buf+ibuf,conf[iz][remap[mu]],sizeof(su3));
		      ibuf+=sizeof(su3);
		    }
		}
      if(ibuf!=sizeof(quad_su3)*glb_vol) crash("did not arrive to the end: %d %d",ibuf,sizeof(quad_su3)*glb_vol);
      
      int nw=fwrite(buf,sizeof(quad_su3),glb_vol,fout);
      if(nw!=glb_vol) crash("did not success in writing");
      nissa_free(buf);
  }
  
  fclose(fout);
  
  nissa_free(conf);
}

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa();
  
  if(rank_tot>1) crash("Cannot run in parallel");
  
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
