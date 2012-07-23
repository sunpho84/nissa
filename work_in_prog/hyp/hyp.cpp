#include "nissa.h"

//level 2 decoration
void hyp_single_rank_dec2_link(su3 link,quad_su3 *conf,int A,int mu,int nu,int rho,double alpha2)
{
  if(nissa_nranks>1) crash("Works only on single rank!");
  if(mu==nu||nu==rho||mu==rho) crash("The three dirs must be different!");
  
  //find the remaining direction
  int eta=0;
  while(eta==mu || eta==nu || eta==rho) eta++;
  
  //take original link
  su3 temp0;
  su3_prod_double(temp0,conf[A][mu],1-alpha2);
  
  //staple and temporary links
  su3 stap,temp1,temp2;
  
  //staple in the positive dir
  int B=loclx_neighup[A][eta];
  int F=loclx_neighup[A][mu];
  unsafe_su3_prod_su3(temp1,conf[A][eta],conf[B][mu]);
  unsafe_su3_prod_su3_dag(stap,temp1,conf[F][eta]);
  
  //staple in the negative dir
  int D=loclx_neighdw[A][eta];
  int E=loclx_neighup[D][mu];
  unsafe_su3_dag_prod_su3(temp1,conf[D][eta],conf[D][mu]);
  unsafe_su3_prod_su3(temp2,temp1,conf[E][eta]);
  su3_summ(stap,stap,temp2);
  
  //summ the two staples with appropriate coef
  su3_summ_the_prod_double(temp0,stap,alpha2/2);
  
  //project the resulting link onto su3
  su3_unitarize_maximal_trace_projecting(link,temp0);
}

//level 1 decoration
void hyp_single_rank_dec1_link(su3 link,quad_su3 *conf,int A,int mu,int nu,double alpha1,double alpha2)
{
  if(nissa_nranks>1) crash("Works only on single rank!");
  if(mu==nu) crash("The two dirs must be different!");
  
  //take original link
  su3 temp0;
  su3_prod_double(temp0,conf[A][mu],1-alpha1);
  
  //reset the staple
  su3 stap,temp1,temp2;
  memset(stap,0,sizeof(su3));
  
  //loop over the third direction
  for(int rho=0;rho<4;rho++)
    if(rho!=mu && rho!=nu)
      {
	//links composing staples
	su3 link1,link2,link3;
	
	//staple in the positive dir
	int B=loclx_neighup[A][rho];
	int F=loclx_neighup[A][mu];
	//prepare the three links
	hyp_single_rank_dec2_link(link1,conf,A,rho,nu,mu,alpha2);
	hyp_single_rank_dec2_link(link2,conf,B,mu,rho,nu,alpha2);
	hyp_single_rank_dec2_link(link3,conf,F,rho,nu,mu,alpha2);
	//compute positive staple
	unsafe_su3_prod_su3(temp1,link1,link2);
	unsafe_su3_prod_su3_dag(temp2,temp1,link3);
	su3_summ(stap,stap,temp2);
	
	//staple in the negative dir
	int D=loclx_neighdw[A][rho];
	int E=loclx_neighup[D][mu];
	//prepare the three links
	hyp_single_rank_dec2_link(link1,conf,D,rho,nu,mu,alpha2);
	hyp_single_rank_dec2_link(link2,conf,D,mu,rho,nu,alpha2);
	hyp_single_rank_dec2_link(link3,conf,E,rho,nu,mu,alpha2);
	//compute negative staple
	unsafe_su3_dag_prod_su3(temp1,link1,link2);
	unsafe_su3_prod_su3(temp2,temp1,link3);
	su3_summ(stap,stap,temp2);
	
	//summ the two staples with appropriate coef
	su3_summ_the_prod_double(temp0,stap,alpha1/4);
      }
  
  //project the resulting link onto su3
  su3_unitarize_maximal_trace_projecting(link,temp0);
}

//level 0 decoration
void hyp_single_rank_dec0_link(su3 link,quad_su3 *conf,int A,int mu,double alpha0,double alpha1,double alpha2)
{
  if(nissa_nranks>1) crash("Works only on single rank!");
  
  //take original link
  su3 temp0;
  su3_prod_double(temp0,conf[A][mu],1-alpha0);
  
  //reset the staple
  su3 stap,temp1,temp2;
  memset(stap,0,sizeof(su3));
  
  //loop over the second direction
  for(int nu=0;nu<4;nu++)
    if(mu!=nu)
      {
	//links composing staples
	su3 link1,link2,link3;
	
	//staple in the positive dir
	int B=loclx_neighup[A][nu];
	int F=loclx_neighup[A][mu];
	//prepare the three links
	hyp_single_rank_dec1_link(link1,conf,A,nu,mu,alpha1,alpha2);
	hyp_single_rank_dec1_link(link2,conf,B,mu,nu,alpha1,alpha2);
	hyp_single_rank_dec1_link(link3,conf,F,nu,mu,alpha1,alpha2);
	//compute positive staple
	unsafe_su3_prod_su3(temp1,link1,link2);
	unsafe_su3_prod_su3_dag(temp2,temp1,link3);
	su3_summ(stap,stap,temp2);
	
	//staple in the negative dir
	int D=loclx_neighdw[A][nu];
	int E=loclx_neighup[D][mu];
	//prepare the three links
	hyp_single_rank_dec1_link(link1,conf,D,nu,mu,alpha1,alpha2);
	hyp_single_rank_dec1_link(link2,conf,D,mu,nu,alpha1,alpha2);
	hyp_single_rank_dec1_link(link3,conf,E,nu,mu,alpha1,alpha2);
	//compute negative staple
	unsafe_su3_dag_prod_su3(temp1,link1,link2);
	unsafe_su3_prod_su3(temp2,temp1,link3);
	su3_summ(stap,stap,temp2);
	
	//summ the two staples with appropriate coef
	su3_summ_the_prod_double(temp0,stap,alpha0/6);
      }
  
  //project the resulting link onto su3
  su3_unitarize_maximal_trace_projecting(link,temp0);
}

//smear a conf on a single rank
void hyp_single_rank_smear_conf(quad_su3 *sm_conf,quad_su3 *conf,double alpha0,double alpha1,double alpha2)
{
  if(nissa_nranks>1) crash("Works only on single rank!");
  
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int mu=0;mu<4;mu++)
      hyp_single_rank_dec0_link(sm_conf[ivol][mu],conf,ivol,mu,alpha0,alpha1,alpha2);
}

int main(int narg,char **arg)
{
  init_nissa();
  
  //set the lattice grid
  int T=8;
  int L=4;
  
  //init the grid
  init_grid(T,L);
  
  //read conf
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  read_ildg_gauge_conf(conf,"../../data/L4T8conf");
  communicate_lx_quad_su3_borders(conf);
  
  //smearing parameters
  double alpha0=0.75;
  double alpha1=0.6;
  double alpha2=0.3;
  
  //smear the conf with hyp
  double hyp_time=-take_time();
  quad_su3 *hyp_smeared_conf=nissa_malloc("hyp_smeared_conf",loc_vol+bord_vol,quad_su3);  
  if(nissa_nranks==1) hyp_single_rank_smear_conf(hyp_smeared_conf,conf,alpha0,alpha1,alpha2);
  hyp_time+=take_time();
  
  //smear the conf with hyp using fast routine
  double hypf_time=-take_time();
  quad_su3 *hyp_fast_smeared_conf=nissa_malloc("hyp_fast_smeared_conf",loc_vol+bord_vol,quad_su3);  
  hyp_smear_conf(hyp_fast_smeared_conf,conf,alpha0,alpha1,alpha2);
  hypf_time+=take_time();
  
  //smeare the conf with APE
  quad_su3 *ape_smeared_conf=nissa_malloc("ape_smeared_conf",loc_vol+bord_vol,quad_su3);  
  ape_spatial_smear_conf(ape_smeared_conf,conf,0.7,1);
  
  //compute plaquette
  master_printf("orig: plaq: %.18g, var: %.18g\n",global_plaquette_lx_conf(conf),global_plaquette_variance_lx_conf(conf));
  if(nissa_nranks==1) master_printf("hyp:  plaq: %.18g, var: %.18g\n",global_plaquette_lx_conf(hyp_smeared_conf),global_plaquette_variance_lx_conf(hyp_smeared_conf));
  else            master_printf("hyp:  plaq: %.18g, var: %.18g\n",0.930002530940892247,0.0396775676154377671);
  master_printf("hypf: plaq: %.18g, var: %.18g\n",global_plaquette_lx_conf(hyp_fast_smeared_conf),global_plaquette_variance_lx_conf(hyp_fast_smeared_conf));
  master_printf("ape:  plaq: %.18g, var: %.18g\n",global_plaquette_lx_conf(ape_smeared_conf),global_plaquette_variance_lx_conf(ape_smeared_conf));
  
  master_printf("Timings: ordinary=%lg s, fast=%lg s.\n",hyp_time,hypf_time);
  
  //free
  nissa_free(conf);
  nissa_free(ape_smeared_conf);
  nissa_free(hyp_smeared_conf);
  nissa_free(hyp_fast_smeared_conf);
  
  close_nissa();
  
  return 0;
}
