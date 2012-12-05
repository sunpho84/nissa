#include "nissa.h"
#include <math.h>

int L=24,T=48;

int main(int narg,char **arg)
{
  init_nissa();
  if(narg<3) crash("Use: %s file_out file_in",arg[0]);
  
  init_grid(T,L);
  
  //open the file and allocate the conf
  FILE *fin=open_file(arg[2],"r");
  quad_su3 *conf=nissa_malloc("conf",loc_vol,quad_su3);
  
  //seek after 4 bytes
  int rc=fseek(fin,4,SEEK_SET);
  if(rc!=0) crash("error while seeking, received: %d",rc);
  
  //read V quad_su3
  rc=fread((void*)conf,sizeof(quad_su3),loc_vol,fin);
  if(rc!=loc_vol) crash("error while reading, received: %d",rc);
  
  //read last entry, that is, the size
  int last;
  rc=fread((void*)&last,sizeof(int),1,fin);
  if(rc!=1) crash("error while reading last byte, obtained %d",rc);
  if(last!=loc_vol*sizeof(quad_su3)) crash("error while reading the conf");

  //transpose the colors
  nissa_loc_vol_loop(ivol)
    for(int mu=0;mu<4;mu++)
      {
	//check su3ity
	complex c;
	su3_det(c,conf[0][0]);
	c[0]-=1;
	double r=sqrt(c[0]*c[0]+c[1]*c[1]);
	if(r>=1.e-15) crash("link %d,%d not in su3, determinant: %lg",ivol,mu,r);
      }
  
  //reorder the conf and write the plaquette
  reorder_read_ildg_gauge_conf(conf);
  master_printf("Plaquette: %lg\n",global_plaquette_lx_conf(conf));
  
  //write the converted conf
  write_ildg_gauge_conf(arg[1],conf,64);

  nissa_free(conf);
  
  fclose(fin);
  
  close_nissa();
  
  return 0;
}
