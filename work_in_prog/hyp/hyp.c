#include "nissa.h"

/*
typedef quad_su3 dec1_quad_su3[3];
typedef dec1_quad_su3 dec2_quad_su3[2];
*/

//level 2 decoration
void hyp_single_rank_dec2_link(su3 link,quad_su3 *conf,int A,int mu,int nu,int rho,double alpha2)
{
  if(rank_tot>1) crash("Works only on single rank!");
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
  if(rank_tot>1) crash("Works only on single rank!");
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
  if(rank_tot>1) crash("Works only on single rank!");
  
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
  if(rank_tot>1) crash("Works only on single rank!");
  
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int mu=0;mu<4;mu++)
      hyp_single_rank_dec0_link(sm_conf[ivol][mu],conf,ivol,mu,alpha0,alpha1,alpha2);
}

//smear a conf on multiple ranks
void hyp_many_rank_smear_conf(quad_su3 *sm_conf,quad_su3 *conf,double alpha0,double alpha1,double alpha2)
{
  //fill the dec2 remapping table
  int dec2_remap_index[4][4][4];
  int idec2_remap=0;
  for(int mu=0;mu<4;mu++)
    for(int nu=0;nu<4;nu++)
      for(int rho=0;rho<4;rho++)
	  if(mu==nu|| mu==rho || nu==rho) dec2_remap_index[mu][nu][rho]=-1;
	  else                            dec2_remap_index[mu][nu][rho]=idec2_remap++;
  
  //fill the dec1 remapping table
  int dec1_remap_index[4][4];
  int idec1_remap=0;
  for(int mu=0;mu<4;mu++)
    for(int nu=0;nu<4;nu++)
      if(mu==nu) dec1_remap_index[mu][nu]=-1;
      else       dec1_remap_index[mu][nu]=idec1_remap++;
  
  //communicate borders of original conf
  communicate_lx_gauge_borders(conf);
  
  /////////////////////////////////////// second level decoration /////////////////////////////////
  
  //allocate dec2 conf
  su3 *dec2_conf[24];
  for(int idec2=0;idec2<24;idec2++) dec2_conf[idec2]=nissa_malloc("dec2_conf",loc_vol+loc_bord,su3);
  
  //loop over external index
  for(int mu=0;mu<4;mu++)
    //loop over the first decoration index
    for(int nu=0;nu<4;nu++)
      if(nu!=mu)
	//loop over the second decoration index
	for(int rho=0;rho<4;rho++)
	  if(rho!=mu && rho!=nu)
	    {
	      //find the remaining direction
	      int eta=0;
	      while(eta==mu || eta==nu || eta==rho) eta++;
	      
	      //find the remapped index
	      int ire0=dec2_remap_index[mu][nu][rho];
	      
	      //loop over local volume
	      for(int A=0;A<loc_vol;A++)
		{
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
		  su3_unitarize_maximal_trace_projecting(dec2_conf[ire0][A],temp0);
		}
	      
	      //communicate borders for future usage
	      communicate_lx_su3_borders(dec2_conf[ire0]);
	    }

  for(int ivol=0;ivol<loc_vol+loc_bord;ivol++)
    for(int ire0=0;ire0<24;ire0++)
      master_printf("%d  %d %lg\n",ivol,ire0,real_part_of_trace_su3_prod_su3_dag(dec2_conf[ire0][ivol],dec2_conf[ire0][ivol]));
  
  /////////////////////////////////////// first level decoration /////////////////////////////////

  //allocate dec1 conf
  su3 *dec1_conf[12];
  for(int idec1=0;idec1<12;idec1++) dec1_conf[idec1]=nissa_malloc("dec1_conf",loc_vol+loc_bord,su3);
  
  //loop over external index
  for(int mu=0;mu<4;mu++)
    //loop over the first decoration index
    for(int nu=0;nu<4;nu++)
      if(nu!=mu)
	//loop over local volume
	for(int A=0;A<loc_vol;A++)
	  {
	    //take original link
	    su3 temp0;
	    su3_prod_double(temp0,conf[A][mu],1-alpha1);
	    
	    //reset the staple
	    su3 stap;
	    memset(stap,0,sizeof(su3));
	    
	    //find the remapped index
	    int ire0=dec1_remap_index[mu][nu];
	    
	    //loop over the second decoration index
	    for(int rho=0;rho<4;rho++)
	      if(rho!=mu && rho!=nu)
		{
		  su3 temp1,temp2;
		  
		  //find the two remampped indices
		  int ire1=dec2_remap_index[rho][nu][mu];
		  int ire2=dec2_remap_index[mu][rho][nu];
		  
		  //staple in the positive dir
		  int B=loclx_neighup[A][rho];
		  int F=loclx_neighup[A][mu];
		  unsafe_su3_prod_su3(temp1,dec2_conf[ire1][A],dec2_conf[ire2][B]);
		  unsafe_su3_prod_su3_dag(temp2,temp1,dec2_conf[ire1][F]);
		  su3_summ(stap,stap,temp2);
		  
		  //staple in the negative dir
		  int D=loclx_neighdw[A][rho];
		  int E=loclx_neighup[D][mu];
		  unsafe_su3_dag_prod_su3(temp1,dec2_conf[ire1][D],dec2_conf[ire2][D]);
		  unsafe_su3_prod_su3(temp2,temp1,dec2_conf[ire1][E]);
		  su3_summ(stap,stap,temp2);
		  
		  //summ the two staples with appropriate coef
		  su3_summ_the_prod_double(temp0,stap,alpha1/4);
		  
		  //project the resulting link onto su3
		  su3_unitarize_maximal_trace_projecting(dec1_conf[ire0][A],temp0);
		}
	    
	    //communicate borders for future usage
	    communicate_lx_su3_borders(dec1_conf[ire0]);
	  }
  
  //free dec2
  for(int idec2=0;idec2<24;idec2++) nissa_free(dec2_conf[idec2]);
  
  /////////////////////////////////////// first level decoration /////////////////////////////////
  
  //loop over external index
  for(int mu=0;mu<4;mu++)
    //loop over local volume
    for(int A=0;A<loc_vol;A++)
      {
	//take original link
	su3 temp0;
	su3_prod_double(temp0,conf[A][mu],1-alpha0);
	
	//reset the staple
	su3 stap;
	memset(stap,0,sizeof(su3));
	
	//loop over the first decoration index
	for(int nu=0;nu<4;nu++)
	  if(nu!=mu)
	    {
	      su3 temp1,temp2;
	      
	      //find the two remampped indices
	      int ire1=dec1_remap_index[nu][mu];
	      int ire2=dec1_remap_index[mu][nu];
	      
	      //staple in the positive dir
	      int B=loclx_neighup[A][nu];
	      int F=loclx_neighup[A][mu];
	      unsafe_su3_prod_su3(temp1,dec1_conf[ire1][A],dec1_conf[ire2][B]);
	      unsafe_su3_prod_su3_dag(temp2,temp1,dec1_conf[ire1][F]);
	      su3_summ(stap,stap,temp2);
	      
	      //staple in the negative dir
	      int D=loclx_neighdw[A][nu];
	      int E=loclx_neighup[D][mu];
	      unsafe_su3_dag_prod_su3(temp1,dec1_conf[ire1][D],dec1_conf[ire2][D]);
	      unsafe_su3_prod_su3(temp2,temp1,dec1_conf[ire1][E]);
	      su3_summ(stap,stap,temp2);
	      
	      //summ the two staples with appropriate coef
	      su3_summ_the_prod_double(temp0,stap,alpha0/6);
	      
	      //project the resulting link onto su3
	      su3_unitarize_maximal_trace_projecting(sm_conf[A][mu],temp0);
	    }
      }
  
  //free dec1
  for(int idec1=0;idec1<12;idec1++) nissa_free(dec1_conf[idec1]);
}

  
//wrapper
void hyp_smear_conf(quad_su3 *sm_conf,quad_su3 *conf,double alpha0,double alpha1,double alpha2)
{
  if(rank_tot==1) hyp_single_rank_smear_conf(sm_conf,conf,alpha0,alpha1,alpha2);
    else crash("not yet implemented in parallel");
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
  quad_su3 *conf=nissa_malloc("conf",loc_vol+loc_bord,quad_su3);  
  read_ildg_gauge_conf(conf,"../../data/L4T8conf");
  communicate_lx_gauge_borders(conf);
  
  //smearing parameters
  double alpha0=0.75;
  double alpha1=0.6;
  double alpha2=0.3;
  
  //smear the conf with hyp
  double hyp_time=-take_time();
  quad_su3 *hyp_smeared_conf=nissa_malloc("hyp_smeared_conf",loc_vol+loc_bord,quad_su3);  
  if(rank_tot==1) hyp_smear_conf(hyp_smeared_conf,conf,alpha0,alpha1,alpha2);
  hyp_time+=take_time();
  
  //smear the conf with hyp using fast routine
  double hypf_time=-take_time();
  quad_su3 *hyp_fast_smeared_conf=nissa_malloc("hyp_fast_smeared_conf",loc_vol+loc_bord,quad_su3);  
  hyp_many_rank_smear_conf(hyp_fast_smeared_conf,conf,alpha0,alpha1,alpha2);
  hypf_time+=take_time();
  
  //smeare the conf with APE
  quad_su3 *ape_smeared_conf=nissa_malloc("ape_smeared_conf",loc_vol+loc_bord,quad_su3);  
  ape_smear_conf(ape_smeared_conf,conf,0.7,1);
  
  //compute plaquette
  master_printf("orig: plaq: %.18g, var: %.18g\n",global_plaquette_lx_conf(conf),global_plaquette_variance_lx_conf(conf));
  if(rank_tot==1) master_printf("hyp:  plaq: %.18g, var: %.18g\n",global_plaquette_lx_conf(hyp_smeared_conf),global_plaquette_variance_lx_conf(hyp_smeared_conf));
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
