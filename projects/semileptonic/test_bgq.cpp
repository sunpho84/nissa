#include "nissa.h"
#include <math.h>
#include "../../src/bgq/hopping_matrix_bgq.h"

void hopping_matrix_expand_to_Q_and_summ_diag_term_bgq_binded(bi_spincolor *out,double kappa,double mu,bi_spincolor *in);
void bgq_Wilson_hopping_matrix_T_VN_comm_and_buff_fill();

const int nbench=100,nbench_port=10;

int seed=100;
int L=24,T=L*2;
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
  
  GET_THREAD_ID();
  
  //timing mpi timing
  double ave_timing=0,var_timing=0;
  for(int ibench=0;ibench<nbench;ibench++)
    {
      double iter_timing=-take_time();
      iter_timing+=take_time();
      
      ave_timing+=iter_timing;
      var_timing+=iter_timing*iter_timing;     
    }
  ave_timing/=nbench;
  var_timing/=nbench;
  var_timing-=ave_timing*ave_timing;
  var_timing=sqrt(var_timing);
  master_printf("timing overherad: %lg +- %lg sec per timing\n",ave_timing,var_timing);
  
  //benchmark thread_barrier
  int nbench_barr=1;
  double barr_time=-take_time();
  for(int ibench=0;ibench<nbench_barr;ibench++)
    THREAD_BARRIER();
  barr_time+=take_time();
  barr_time/=nbench_barr;
  //sleep(1);
  //and pure omp barrier
  double omp_barr_time=-take_time();
  if(0)
    for(int ibench=0;ibench<nbench_barr;ibench++)
      {
#pragma omp barrier
	{
	  int i=1;
	}
	master_printf("barr %d passed\n",ibench);
      }
  omp_barr_time+=take_time();
  omp_barr_time/=nbench_barr;
  master_printf("thread barrier: %lg sec/barr, %lg sec/omp_barr\n",barr_time,omp_barr_time);

  //benchmark communications in the 8 dirs
  for(int jdir=0;jdir<8;jdir++)
    {
      //set dir
      int comm_dir[8];
      for(int kdir=0;kdir<8;kdir++) comm_dir[kdir]=0;
      comm_dir[jdir]=1;
      
      int idir=jdir%4;
      if(paral_dir[idir])
	{
	  //take time of n comms
	  double time=-take_time();
	  for(int ibench=0;ibench<nbench;ibench++)
	    {
	      comm_start(lx_spincolor_comm,comm_dir,sizeof(spincolor)*bord_dir_vol[idir]);
	      comm_wait(lx_spincolor_comm);
	    }
	  time+=take_time();
	  
	  //print out
	  double size=bord_dir_vol[idir]*sizeof(spincolor)/1024.0/1024*nbench;
	  if(rank==0 && thread_id==0) printf("Rank %d, speed in %d dir: %lg Mb / %lg sec = %lg Mb/sec\n",rank,idir,size,time,size/time);
	}
    }
  
  //apply a fixed number of time
  spincolor *out=nissa_malloc("out",loc_vol+bord_vol,spincolor);
  spincolor *tmp=nissa_malloc("tmp",loc_vol+bord_vol,spincolor);
  double port_time=-take_time();
  for(int ibench=0;ibench<nbench_port;ibench++)
    apply_tmQ2(out,conf,kappa,tmp,mu,in);
  port_time+=take_time();
  port_time/=nbench_port;
  
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
  for(int ibench=0;ibench<nbench;ibench++)
    apply_tmQ2_bgq(bi_out,bi_conf,kappa,mu,bi_in);
  bgq_time+=take_time();
  bgq_time/=nbench;
  
  //unmap to compare
  spincolor *un_out=nissa_malloc("un_out",loc_vol+bord_vol,spincolor);
  bgqlx_spincolor_remap_to_lx(un_out,bi_out);
  
  //compute average diff
  double diff;
  double_vector_subt((double*)un_out,(double*)un_out,(double*)out,loc_vol*sizeof(spincolor)/sizeof(double));
  double_vector_glb_scalar_prod(&diff,(double*)un_out,(double*)un_out,loc_vol*sizeof(spincolor)/sizeof(double));
  master_printf("Diff: %lg\n",diff);

  //benchmark pure hopping matrix application
  double hop_bgq_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    apply_Wilson_hopping_matrix_bgq_binded_nocomm_nobarrier(bi_conf,0,loc_volh,bi_in);
  hop_bgq_time+=take_time();
  hop_bgq_time/=nbench;
  int nflops_hop=1152*loc_vol;
  master_printf("hop_bgq_time: %lg sec, %d flops, %lg Mflops\n",hop_bgq_time,nflops_hop,nflops_hop*1e-6/hop_bgq_time);
  
  //benchmark expansion
  double exp_bgq_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    hopping_matrix_expand_to_Q_and_summ_diag_term_bgq_binded(bi_out,kappa,mu,bi_in);
  exp_bgq_time+=take_time();
  exp_bgq_time/=nbench;
  int nflops_exp=312*loc_vol;
  master_printf("exp_bgq_time: %lg sec, %d flops, %lg Mflops\n",exp_bgq_time,nflops_exp,nflops_exp*1e-6/exp_bgq_time);
  master_printf("(hop+exp)_bgq_time: %lg sec, %d flops, %lg Mflops\n",
		hop_bgq_time+exp_bgq_time,nflops_exp+nflops_hop,(nflops_exp+nflops_hop)*
		1e-6/(hop_bgq_time+exp_bgq_time));
  
  //benchmark pure spi communication without data moving
  double pure_spi_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    {
      comm_start(lx_halfspincolor_comm);
      comm_wait(lx_halfspincolor_comm);
    }
  pure_spi_time+=take_time();
  pure_spi_time/=nbench;
  double data_exch=lx_halfspincolor_comm.tot_mess_size/1024.0/1024.0;
  master_printf("pure_comm_time: %lg sec, %lg Mbytes, %lg Mb/sec\n",pure_spi_time,data_exch,data_exch/pure_spi_time);
  
  //benchmark buff filling
  double buff_filling_time=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    bgq_Wilson_hopping_matrix_T_VN_comm_and_buff_fill();
  buff_filling_time+=take_time();
  buff_filling_time/=nbench;
  master_printf("buff_filling_time: %lg sec\n",buff_filling_time);
  
  nissa_free(bi_conf);
  nissa_free(bi_in);
  nissa_free(bi_out);
  nissa_free(un_out);
#endif

  ////////////////////////////////
  
  master_printf("In+conf+addr size: %lg Mbytes (L2 cache: 32 Mb)\n",
		(loc_vol*(sizeof(quad_su3)+sizeof(spincolor))+8*loc_volh*sizeof(bi_halfspincolor*))/1024.0/1024.0);
  master_printf("Time to apply %d time:\n",nbench);
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
