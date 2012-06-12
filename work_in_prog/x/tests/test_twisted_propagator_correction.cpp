#include <math.h>

#include "nissa.h"

#include "../src/propagators/twisted_propagator.h"
#include "../src/propagators/tlSym_gluon_propagator.h"
#include "../src/propagators/Wilson_gluon_propagator.h"
#include "../src/types/types_routines.h"
#include "../src/routines/fourier.h"
#include "../src/routines/shift.h"
#include "../src/diagrams/propagator_self_energy.h"
#include "../src/diagrams/tadpole.h"

spinspin *corr2_x,*corr2_p;
spinspin *corr1_x,*corr1_p;

int comp=0*1;
int L;

double pana_c1_id_corr(double a2p2,double eps,double alpha)
{
  double a=eps+1.529131270100*alpha;
  double b=-48.932201;
  return (a2p2*a+b)/(16*M_PI*M_PI);
}

double pana_c2_id_corr(double a2p2,double eps,double alpha)
{
  double a=eps-5.39301570*alpha-0.5*(3-2*alpha)*log(a2p2);
  double b=-2.502511;
  
  return (a2p2*a+b)/(16*M_PI*M_PI);
}

double pana_c2_ipslash_corr(double a2p2,double eps,double alpha)
{return (eps-7.850272109*alpha+alpha*log(a2p2))/(16*M_PI*M_PI);}

double real_part_of_trace_with_sinpslash(spinspin *q,int imom,quark_info qu)
{
  complex tr={0,0};
  double a2pt2=0,a2p2=0;
  for(int mu=0;mu<4;mu++)
    {
      int nu=nissa_map_mu[mu];
      double ap=M_PI*(2*glb_coord_of_loclx[imom][mu]+qu.bc[mu])/glb_size[mu];
      double apt=2*sin(ap/2);
      a2pt2+=apt*apt;
      a2p2+=ap*ap;
      //print_dirac(base_gamma+nu);
      for(int id=0;id<4;id++)
	{
	  complex t;
	  int ig=base_gamma[nu].pos[id];
	  unsafe_complex_prod(t,base_gamma[nu].entr[id],q[imom][ig][id]);
	  complex_summ_the_prod_double(tr,t,ap);
	  //printf("mu=%d nu=%d ig=%d id=%d ap=%lg ",mu,nu,ig,id,ap);
	  //printf("pr=(%lg,%lg) ",q[imom][ig][id][0],q[imom][ig][id][1]);
	  //printf("g=(%lg,%lg) ",base_gamma[nu].entr[id][0],base_gamma[nu].entr[id][1]);
	  //printf("contr=(%lg,%lg) ",t[0],t[1]);
	  //printf("tot=(%lg,%lg)\n",tr[0],tr[1]);
	}
      //master_printf("mu=%d %lg\n",mu,tr[1]);
    }
  return tr[1]/4/a2pt2;
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
void init_test()
{
  //init the grid
  init_grid(L,L);
  
  //allocate corrections
  corr1_x=nissa_malloc("corr1_x",loc_vol,spinspin);
  corr1_p=nissa_malloc("corr1_p",loc_vol,spinspin);
  corr2_x=nissa_malloc("corr2_x",loc_vol,spinspin);
  corr2_p=nissa_malloc("corr2_p",loc_vol,spinspin);
}

//close the program
void close_test()
{
  nissa_free(corr2_x);
  nissa_free(corr2_p);
  nissa_free(corr1_p);
  nissa_free(corr1_x);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa();
  
  if(narg<2) crash("use: %s L",arg[0]);
  L=atoi(arg[1]);
  
  init_test();
  
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
  
  compute_tadpole_diagram_in_mom_space(corr1_p,qu,gl);
  compute_tadpole_diagram_in_x_space(corr1_x,gl);
  pass_spinspin_from_x_to_mom_space(corr1_x,corr1_x,qu.bc);
  
  /////////////////////////////////// correction D2 ///////////////////////////////////
  
  compute_self_energy_twisted_propagator_in_x_space(corr2_x,qu,gl);
  pass_spinspin_from_x_to_mom_space(corr2_x,corr2_x,qu.bc);
  
  ///////////////////////////////// correction in P ////////////////////////////  
  
  if(comp && rank_tot==1)
    {
      compute_self_energy_twisted_propagator_in_mom_space(corr2_p,qu,gl);

      double tt=0;
      nissa_loc_vol_loop(ivol)
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
	    master_printf("%d %lg\n",ivol,t2);
	    print_spinspin(temp);
	    master_printf("x-space: \n");
	    print_spinspin(corr2_x[ivol]);
	    master_printf("p-space: \n");
	    print_spinspin(corr2_p[ivol]);
	  }
      }
      tt=glb_reduce_double(tt);
      master_printf("\nDifference between mom and x-space computation: %lg\n\n",sqrt(tt/glb_vol));
    }
  
  //////////////////////////////////// output ////////////////////////////

  //   coefficients for tad-pole diagram
  double eps1_11=4.5873938103; //plaquette
  //double eps1_11=3.535587351;  //tlSym
  
  //   coefficients for self-energy diagram
  double eps2_01=7.4696262;
  double eps2_11=8.2395316;//plaquette
  //double eps2_11=7.16084231;//tlSym
  
  char path_c1[1024],path_c2[1024];
  sprintf(path_c1,"outc1_%02d_d",L);
  sprintf(path_c2,"outc2_%02d_d",L);
  FILE *fout_c1=open_file(path_c1,"w");
  FILE *fout_c2=open_file(path_c2,"w");
  nissa_loc_vol_loop(imom)
    {
      int w=1;
      double a2pt2=0,a2p2=0,th=0;
      for(int mu=0;mu<4;mu++)
	{
	  double ap=M_PI*(2*glb_coord_of_loclx[imom][mu]+qu.bc[mu])/glb_size[mu];
	  double apt=sin(ap);
	  a2p2+=ap*ap;
	  a2pt2+=apt*apt;
	  w=w&&(glb_coord_of_loclx[imom][mu]<=glb_size[mu]/4);
	  th+=ap/2;
	}
      th=acos(th/sqrt(a2p2));
      
      double c1_tp2_x=real_part_of_trace_with_id(corr1_x,imom,qu)*glb_vol;
      double c1_tp2=real_part_of_trace_with_id(corr1_p,imom,qu)*glb_vol;
      double c2_tp2=real_part_of_trace_with_id(corr2_x,imom,qu)*glb_vol;
      double c2_tpsl=real_part_of_trace_with_sinpslash(corr2_x,imom,qu)*glb_vol;
      
      double pa_c1_id=pana_c1_id_corr(a2p2,eps1_11,alpha);
      double pa_c2_id=pana_c2_id_corr(a2p2,eps2_11,alpha);
      double pa_c2_ips=pana_c2_ipslash_corr(a2pt2,eps2_01,alpha);
      if(fabs(th)<0.2 && w)
	{
	  fprintf(fout_c1,"%lf %lg\t%lf %lg %lg\t%lg\n",a2p2,pa_c1_id,a2pt2,c1_tp2_x,c1_tp2,th);
	  fprintf(fout_c2,"%lf\t%lg\t%lg\t%lg\t%lg\n",a2p2,pa_c2_id,a2pt2,c2_tp2,th);
	}
    }
  fclose(fout_c1);
  fclose(fout_c2);
  
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
  master_printf("a2p2: %lg\n",a2p2);
  master_printf("att: %lg\n",pana_c2_id_corr(a2p2,eps2_11,alpha));
  
  if(rank_tot==1 && comp)
    {
      master_printf("p-space: \n");
      print_spinspin(corr2_p[lx]);
      master_printf("\n");
    }
  
  master_printf("x-space: \n");
  if(rank==rx) print_spinspin(corr2_x[lx]);
  
  close_test();
  
  return 0;
}
