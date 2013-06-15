#include "nissa.h"

int seed=100;
int napp=100;
int L=16,T=L*2;
coords source_coord={0,0,0,0};
double mu=0.03,kappa=0.137;

void in_main(int narg,char **arg)
{
  //init
  init_grid(T,L); 
  start_loc_rnd_gen(seed);
  
  //create random conf
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  nissa_loc_vol_loop(ivol)
    for(int mu=0;mu<4;mu++)
      su3_put_to_rnd(conf[ivol][mu],loc_rnd_gen[ivol]);
  set_borders_invalid(conf);
  
  //create random in vector
  spincolor *in=nissa_malloc("in",loc_vol+bord_vol,spincolor);
  generate_undiluted_source(in,RND_Z4,-1);
  
  //benchmark communications in the 8 dirs
  double time[8];
  for(int idir=0;idir<8;idir++)
    if(paral_dir[idir%4])
      {
	//set dir
	int comm_dir[8];
	for(int jdir=0;jdir<8;jdir++) comm_dir[jdir]=0;
	comm_dir[idir]=1;
	
	//take time of n comms
	int ncomm=100;
	time[idir]=-take_time();
	for(int ibench=0;ibench<ncomm;ibench++)
	  {
	    comm_start(lx_spincolor_comm,comm_dir);
	    comm_wait(lx_spincolor_comm);
	  }
	time[idir]+=take_time();
	
	//print out
	double size=bord_dir_vol[idir%4]*sizeof(spincolor)/1024.0/1024*ncomm;
	printf("Rank %d, speed in %d dir: %lg Mb / %lg sec = %lg Mb/sec\n",rank,idir,size,time[idir%4],size/time[idir]);
      }
  
  //apply a fixed number of time
  spincolor *out=nissa_malloc("out",loc_vol+bord_vol,spincolor);
  spincolor *tmp=nissa_malloc("tmp",loc_vol+bord_vol,spincolor);
  double port_time=-take_time();
  for(int iapp=0;iapp<napp;iapp++)
    apply_tmQ2(out,conf,kappa,tmp,mu,in);
  port_time+=take_time();
  
  //////////////////////////////////////

#ifdef EXP_BGQ
  
  //remap conf to bgq
  bi_oct_su3 *bi_conf=nissa_malloc("bi_conf",loc_vol+bord_vol,bi_oct_su3);
  lx_conf_remap_to_bgqlx(bi_conf,conf);
  
  //remap in to bgq
  bi_spincolor *bi_in=nissa_malloc("bi_in",loc_vol+bord_vol,bi_spincolor);
  lx_spincolor_remap_to_bgqlx(bi_in,in);
  
  //apply bgq
  bi_spincolor *bi_out=nissa_malloc("bi_out",loc_vol+bord_vol,bi_spincolor);
  double bgq_time=-take_time();
  for(int iapp=0;iapp<napp;iapp++)
    apply_tmQ2_bgq(bi_out,bi_conf,kappa,mu,bi_in);
  bgq_time+=take_time();
  
  //unmap to compare
  spincolor *un_out=nissa_malloc("un_out",loc_vol+bord_vol,spincolor);
  bgqlx_spincolor_remap_to_lx(un_out,bi_out);
  
  //compute average diff
  double diff;
  double_vector_subt((double*)un_out,(double*)un_out,(double*)out,loc_vol*sizeof(spincolor)/sizeof(double));
  double_vector_glb_scalar_prod(&diff,(double*)un_out,(double*)un_out,loc_vol*sizeof(spincolor)/sizeof(double));
  master_printf("Diff: %lg\n",diff);
  
  nissa_free(bi_conf);
  nissa_free(bi_in);
  nissa_free(bi_out);
  nissa_free(un_out);
  
#endif

  ////////////////////////////////
  
  master_printf("Time to apply %d time:\n",napp);
  master_printf(" %lg sec in port mode\n",port_time);
#ifdef EXP_BGQ
  master_printf(" %lg sec in bgq mode\n",bgq_time);
#endif

  nissa_free(conf);
  nissa_free(in);
  nissa_free(tmp);
  nissa_free(out);
  
  close_nissa();
}
 
int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  
  return 0;
}
