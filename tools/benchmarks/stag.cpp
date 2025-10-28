#include <nissa.hpp>

using namespace nissa;

void incr(double &sx,double &s2x,double init,double fina)
{
  double used=fina-init;
  sx+=used;
  s2x+=used*used;
}

void analyse(double &sx,double &s2x,int nbench)
{
  sx/=nbench;
  s2x/=nbench;
  s2x-=sx*sx;
  s2x=sqrt(s2x);
}

//print the statistic
void print_stat(const char *what,double time,int n,int flops)
{
  MASTER_PRINTF("time to %s %d times: %lg s, %lg per iter",what,n,time,time/std::max(n,1));
  MASTER_PRINTF(", %lg MFlop/s\n",flops*1e-6*n/(time?time:1));
}

void bench_stag()
{
  CRASH(""); 
  // //conf
  // quad_su3 *conf[2];
  // for(int eo=0;eo<2;eo++) conf[eo]=nissa_malloc("conf",loc_volh+bord_volh,quad_su3);
  // generate_hot_eo_conf(conf);
  // oct_su3 *Lebconf[2];
  // for(int eo=0;eo<2;eo++)
  //   {
  //     Lebconf[eo]=nissa_malloc("Lebconf",loc_volh+bord_volh,oct_su3);
  //     remap_loceo_conf_to_Lebeo_oct(Lebconf[eo],conf,eo);
  //   }
  
  // //in
  // color *in=nissa_malloc("in",loc_volh+bord_volh,color);
  // generate_fully_undiluted_eo_source(in, RND_GAUSS,-1,EVN);
  // color *Lebin=nissa_malloc("Lebin",loc_volh+bord_volh,color);
  // remap_loc_ev_or_od_to_Leb_vector(Lebin,in,EVN);
  
  // //temp and out
  // color *temp=nissa_malloc("temp",loc_volh+bord_volh,color);
  // color *out=nissa_malloc("out",loc_volh+bord_volh,color);
  // color *Lebout=nissa_malloc("Lebout",loc_volh+bord_volh,color);
  // color *outrec=nissa_malloc("outrec",loc_volh+bord_volh,color);
  
  // int nbench=500;
  // double mass2=1;
  
  // double t;
  
  // RESET_TIMING(tot_comm_time,ntot_comm);
  
  // MASTER_PRINTF("\n");
  // t=-take_time();
  // for(int ibench=0;ibench<nbench;ibench++)
  //   {
  //     set_borders_invalid(in);
  //     apply_stD2ee_m2(out,conf,temp,mass2,in);
  //   }
  // t+=take_time();
  // print_stat("apply staggered operator",t,nbench,1158*loc_volh);
  
  // t=-take_time();
  // for(int ibench=0;ibench<nbench;ibench++)
  //   {
  //     set_borders_invalid(Lebin);
  //     apply_stD2Leb_ee_m2(Lebout,Lebconf,temp,mass2,Lebin);
  //   }
  // t+=take_time();
  
  // print_stat("apply Leb staggered operator",t,nbench,1158*loc_volh);
  // MASTER_PRINTF("\n");
  
  // remap_Leb_ev_or_od_to_loc_vector(outrec,Lebout,EVN);
  // double_vector_subtassign((double*)outrec,(double*)out,loc_volh*sizeof(color)/sizeof(double));
  // double n2diff=double_vector_glb_norm2(outrec,loc_volh);
  // double n2=double_vector_glb_norm2(out,loc_volh);
  // MASTER_PRINTF("Rel norm of the diff: %lg\n",sqrt(n2diff/n2));
  
  // MASTER_PRINTF("Timing to do %d communications: %lg s, %lg each\n",ntot_comm,tot_comm_time,tot_comm_time/ntot_comm);
  
  // //free
  // for(int eo=0;eo<2;eo++)
  //   {
  //     nissa_free(conf[eo]);
  //     nissa_free(Lebconf[eo]);
  //   }
  // nissa_free(in);
  // nissa_free(Lebin);
  // nissa_free(temp);
  // nissa_free(out);
  // nissa_free(Lebout);
  // nissa_free(outrec);
}

void bench_Wils()
{
  
  //conf
  quad_su3 *conf=nissa_malloc("conf",locVol+bord_vol,quad_su3);
  generate_hot_lx_conf(conf);
  
  //in
  spincolor *in=nissa_malloc("in",locVol+bord_vol,spincolor);
  CRASH("reimplement");
  //generate_fully_undiluted_eo_source(in, RND_GAUSS,-1,EVN);
  
  //temp and out
  spincolor *temp=nissa_malloc("temp",locVol+bord_vol,spincolor);
  spincolor *out=nissa_malloc("out",locVol+bord_vol,spincolor);
  spincolor *outrec=nissa_malloc("outrec",locVol+bord_vol,spincolor);
  
  int nbench=500;
  double mu=1,kappa=0.2;
  
  double t;
  
  RESET_TIMING(tot_comm_time,ntot_comm);
  
  MASTER_PRINTF("\n");
  t=-take_time();
  for(int ibench=0;ibench<nbench;ibench++)
    {
      set_borders_invalid(in);
      apply_tmQ2(out,conf,kappa,temp,mu,in);
    }
  t+=take_time();
  print_stat("apply tm operator",t,nbench,1158*2*locVol);
  
  MASTER_PRINTF("Timing to do %d communications: %lg s, %lg each\n",ntot_comm,tot_comm_time,tot_comm_time/ntot_comm);
  
  //free
  nissa_free(conf);
  nissa_free(in);
  nissa_free(temp);
  nissa_free(out);
  nissa_free(outrec);
}

void in_main(int narg,char **arg)
{
  if(narg<3) CRASH("Use %s L T",arg[0]);
  
  //grid
  int L=atoi(arg[1]),T=atoi(arg[2]);
  initGrid(T,L);
  start_loc_rnd_gen(1234);
  
  bench_stag();
  bench_Wils();
}

int main(int narg,char **arg)
{
  initNissa_threaded(narg,arg,in_main);
  closeNissa();
  
  return 0;
}
