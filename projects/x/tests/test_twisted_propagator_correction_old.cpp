#include <math.h>
#include <stdlib.h>

#include "nissa.hpp"
using namespace std;

#include "../src/propagators/twisted_propagator.hpp"
#include "../src/propagators/tlSym_gluon_propagator.hpp"
#include "../src/propagators/Wilson_gluon_propagator.hpp"
#include "../src/types/types_routines.hpp"
#include "../src/routines/fourier.hpp"
#include "../src/routines/shift.hpp"
#include "../src/diagrams/propagator_self_energy.hpp"
#include "../src/diagrams/tadpole.hpp"

spinspin *corr2_x,*corr2_p;
spinspin *corr1_x,*corr1_p;

int comp_p=0;
int comp_tad=0;
int L;

double pana_c1_id_corr(double a2p2,gluon_info gl)
{
  double alpha=gl.alpha;
  double eps1_11;
  double c1=gl.c1;

  if(c1!=0) eps1_11=3.535587351;  //tlSym
  else eps1_11=4.5873938103; //plaquette

  double a=eps1_11+1.529131270100*alpha;
  double b=-48.932201;
  
  return (a2p2*a+b)/(16*M_PI*M_PI);
}

double pana_c1_ig_corr(double ap_mu,double a2p2,gluon_info gl)
{
  double a2p2_mu=ap_mu*ap_mu;
  double alpha=gl.alpha;
  double eps1_01,eps1_21;
  double c1=gl.c1;
  
  if(c1!=0) //tlSym
    {
      eps1_01=7.071174701;
      eps1_21=-1.1785291169;
    }
  else
    {
      eps1_01=9.174787621;
      eps1_21=-1.5291312701;
    }
      
  double A=eps1_01+3.050262540200*alpha;
  double B=a2p2_mu*(eps1_21-0.509710423367*alpha);
  
  return ap_mu*(A+B)/(16*M_PI*M_PI);
}

double pana_c2_id_corr(double a2p2,gluon_info gl)
{
  double alpha=gl.alpha;
  double eps2_11;
  double c1=gl.c1;
  
  if(c1!=0) eps2_11=7.16084231; //tlSym
  else      eps2_11=8.2395316;  //plaquette
  
  double a=eps2_11-5.39301570*alpha-0.5*(3-2*alpha)*log(a2p2);
  double b=-2.502511;
  
  return (a2p2*a+b)/(16*M_PI*M_PI);
}

double pana_c2_ig_corr(double ap_mu,double a2p2,gluon_info gl)
{
  double a2p2_mu=ap_mu*ap_mu;
  double a4p4_mu=a2p2_mu*a2p2_mu;
  double alpha=gl.alpha;
  double eps2_01,eps2_21,eps_24;
  double c1=gl.c1,C2=c1;
  
  if(c1!=0) //tlSym
    {
      eps2_01=5.95209802;
      eps2_21=-3.0693492;
      eps_24=-1.14716212;
    }
  else
    {
      eps2_01=7.4696262;
      eps2_21=-3.21623339;
      eps_24=-1.5048070;
    }
      
  double A=eps2_01-7.850272109*alpha+alpha*log(a2p2);
  double B=a2p2_mu*(eps2_21+1.016711991*alpha+(101.0/120-11.0/30*C2-alpha/6)*log(a2p2));
  double C=a2p2*(eps_24+1.51604667*alpha+(59.0/240+c1/2+C2/60-0.25*(1.5*alpha))*log(a2p2));
  double D=a4p4_mu/a2p2*(-3.0/80-C2/10-5.0/48*alpha);
  
  return ap_mu*(A+B+C+D)/(16*M_PI*M_PI);
}

double real_part_of_trace_with_igamma(spinspin *q,int imom,int mu,quark_info qu)
{
  complex tr={0,0};
  if(mu<0||mu>=4) CRASH("mu=%d",mu);
  
  int nu=map_mu[mu];
  for(int id=0;id<4;id++)
    {
      complex t;
      int ig=base_gamma[nu].pos[id];
      unsafe_complex_prod(t,base_gamma[nu].entr[id],q[imom][ig][id]);
      complex_summassign(tr,t);
    }
  
  return tr[1]/4;
}

double real_part_of_trace_with_id(spinspin *q,int imom,quark_info qu)
{
  double a2pt2=0,a2p2=0,t=0;
  for(int mu=0;mu<4;mu++)
    {
      double ap=M_PI*(2*glb_coord_of_loclx[imom][mu]+qu.bc[mu])/glb_size[mu];
      double apt=2*sin(ap/2);
      a2p2+=ap*ap;
      a2pt2+=apt*apt;
      t+=q[imom][mu][mu][0];
    }
  
  return t/4;
}

//initialize the program
void init_test(int narg,char **arg)
{
  //init the grid
  init_grid(L,L);
  
  //allocate corrections
  if(comp_tad)corr1_x=nissa_malloc("corr1_x",loc_vol,spinspin);
  if(comp_p && comp_tad) corr1_p=nissa_malloc("corr1_p",loc_vol,spinspin);
  corr2_x=nissa_malloc("corr2_x",loc_vol,spinspin);
  if(comp_p) corr2_p=nissa_malloc("corr2_p",loc_vol,spinspin);
}

//close the program
void close_test()
{
  nissa_free(corr2_x);
  if(comp_p) nissa_free(corr2_p);
  if(comp_tad) nissa_free(corr1_x);
  if(comp_p && comp_tad) nissa_free(corr1_p);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  if(narg<2) CRASH("use: %s L",arg[0]);
  L=atoi(arg[1]);
  
  init_test(narg,arg);
  
  double null_theta[4]={0,0,0,0};
  double small=1.e-6;
  double small_theta[4]={small,small,small,small};
  double rand_theta[4]={0.1,0.3,0.6,0.4};
  double anti_theta[4]={1,0,0,0};
  
  //quark
  double quark_theta[4];memcpy(quark_theta,null_theta,sizeof(double)*4);
  double kappa=1.0/8;
  double mass=0;
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  
  //gluon
  double gluon_theta[4];memcpy(gluon_theta,null_theta,sizeof(double)*4);
  double alpha=1;
  gluon_info gl=create_Wilson_gluon_info(alpha,gluon_theta);
  
  /////////////////////////////////// correction D1 ///////////////////////////////////
  
  if(comp_tad)
    {
      if(comp_p) compute_tadpole_diagram_in_mom_space(corr1_p,qu,gl);
      compute_tadpole_diagram_in_x_space(corr1_x,qu,gl);
      pass_spinspin_from_x_to_mom_space(corr1_x,corr1_x,qu.bc);
    }

  /////////////////////////////////// correction D2 ///////////////////////////////////
  
  compute_self_energy_twisted_diagram_in_x_space(corr2_x,qu,gl);
  pass_spinspin_from_x_to_mom_space(corr2_x,corr2_x,qu.bc);
  
  ///////////////////////////////// correction in P ////////////////////////////  
  
  if(comp_p && nranks==1)
    {
      compute_self_energy_twisted_diagram_in_mom_space(corr2_p,qu,gl);

      double tt=0;
      NISSA_LOC_VOL_LOOP(ivol)
      {
	spinspin temp;
	coords opp;
	for(int mu=0;mu<4;mu++) opp[mu]=(glb_size[mu]-glb_coord_of_loclx[ivol][mu])%glb_size[mu];
	int iopp=glblx_of_coord(opp);
	spinspin_subt(temp,corr2_x[ivol],corr2_p[ivol]);
	double t2=real_part_of_trace_spinspin_prod_spinspin_dag(temp,temp);
	tt+=t2;
	if(fabs(t2)>1.e-10) 
	  {
	    MASTER_PRINTF("%d %lg\n",ivol,t2);
	    spinspin_print(temp);
	    MASTER_PRINTF("x-space: \n");
	    spinspin_print(corr2_x[ivol]);
	    MASTER_PRINTF("p-space: \n");
	    spinspin_print(corr2_p[ivol]);
	  }
      }
      tt=glb_reduce_double(tt);
      MASTER_PRINTF("\nDifference between mom and x-space computation: %lg\n\n",sqrt(tt/glb_vol));
    }
  
  //////////////////////////////////// output ////////////////////////////

  int gtr=1;
  char path_c1_id[1024],path_c2_id[1024];
  char path_c1_ipsl[1024],path_c2_ipsl[1024];
  sprintf(path_c1_id,"outc1_id_%02d_d",L);
  sprintf(path_c1_ipsl,"outc1_ig%d_%02d_d",gtr,L);
  sprintf(path_c2_id,"outc2_id_%02d_d",L);
  sprintf(path_c2_ipsl,"outc2_ig%d_%02d_d",gtr,L);
  FILE *fout_c1_id=open_file(path_c1_id,"w");
  FILE *fout_c1_ipsl=open_file(path_c1_ipsl,"w");
  FILE *fout_c2_id=open_file(path_c2_id,"w");
  FILE *fout_c2_ipsl=open_file(path_c2_ipsl,"w");
  NISSA_LOC_VOL_LOOP(imom)
    {
      int w=1;
      double a2pt2=0,a2p2=0,a4p4=0,th=0;
      double ap[4];
      for(int mu=0;mu<4;mu++)
	{
	  ap[mu]=M_PI*(2*glb_coord_of_loclx[imom][mu]+qu.bc[mu])/glb_size[mu];
	  double apt=sin(ap[mu]);
	  a2p2+=ap[mu]*ap[mu];
	  a4p4+=pow(ap[mu],4);
	  a2pt2+=apt*apt;
	  w=w&&(glb_coord_of_loclx[imom][mu]<=glb_size[mu]/4);
	  th+=ap[mu]/2;
	}
      th=acos(th/sqrt(a2p2));
      
      double c1_tid=(comp_tad)?real_part_of_trace_with_id(corr1_x,imom,qu)*glb_vol:0;
      double c1_tpsl=(comp_tad)?real_part_of_trace_with_igamma(corr1_x,imom,gtr,qu)*glb_vol:0;
      double c2_tid=real_part_of_trace_with_id(corr2_x,imom,qu)*glb_vol;
      double c2_tpsl=real_part_of_trace_with_igamma(corr2_x,imom,gtr,qu)*glb_vol;
      
      double pa_c1_id=pana_c1_id_corr(a2p2,gl);
      double pa_c1_ipsl=pana_c1_ig_corr(ap[gtr],a2p2,gl);
      double pa_c2_id=pana_c2_id_corr(a2p2,gl);
      double pa_c2_ipsl=pana_c2_ig_corr(ap[gtr],a2p2,gl);
      if(fabs(th)<0.2 && w)
	{
	  fprintf(fout_c1_id,  "%lf %lg\t%lf %lg\t%lg\n",a2p2,pa_c1_id,  a2pt2,c1_tid, th);
	  fprintf(fout_c1_ipsl,"%lf %lg\t%lf %lg\t%lg\n",a2p2,pa_c1_ipsl,a2pt2,c1_tpsl,th);
	  fprintf(fout_c2_id,  "%lf %lg\t%lf %lg\t%lg\n",a2p2,pa_c2_id,  a2pt2,c2_tid, th);
	  fprintf(fout_c2_ipsl,"%lf %lg\t%lf %lg\t%lg\n",a2p2,pa_c2_ipsl,a2pt2,c2_tpsl,th);
	}
    }
  fclose(fout_c1_id);
  fclose(fout_c1_ipsl);
  fclose(fout_c2_id);
  fclose(fout_c2_ipsl);
  
  //compute the point to be printed
  coords ix={0,0,0,1};
  int lx,rx;
  get_loclx_and_rank_of_coord(&lx,&rx,ix);
  
  double a2p2=0;
  for(int mu=0;mu<4;mu++)
    {
      double p=M_PI*(2*ix[mu]+qu.bc[mu])/glb_size[mu];
      a2p2+=p*p;
    }
  MASTER_PRINTF("a2p2: %lg\n",a2p2);
  MASTER_PRINTF("att: %lg\n",pana_c2_id_corr(a2p2,gl));
  
  if(nranks==1 && comp_p)
    {
      MASTER_PRINTF("p-space: \n");
      spinspin_print(corr2_p[lx]);
      MASTER_PRINTF("\n");
    }
  
  MASTER_PRINTF("x-space: \n");
  if(rank==rx) spinspin_print(corr2_x[lx]);
  
  close_test();
  
  return 0;
}
