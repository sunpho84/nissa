#include "nissa.h"

//level 2 decoration
void hyp_single_rank_dec2_link(su3 link,quad_su3 *conf,int A,int mu,int nu,int rho,double alpha2)
{
  if(nrank>0) crash("Works only on single rank!");
  if(mu==nu||nu==rho||mu==rho) crash("The three dirs must be different!");
  
  //find the remaining direction
  int eta=0;
  while(eta==mu || eta==nu || eta==rho) eta++;
  
  //take original link
  su3 temp0;
  su3_prod_double(temp0,conf[A][mu],1-alpha2);
  
  //reset the staple
  su3 stap,temp1,temp2;
  memset(stap,0,sizeof(su3));
  
  //staple in the positive dir
  int B=loclx_neighup[A][eta];
  int F=loclx_neighup[A][mu];
  unsafe_su3_prod_su3(temp1,conf[A][eta],conf[B][mu]);
  unsafe_su3_prod_su3_dag(temp2,temp1,conf[F][eta]);
  su3_summ(stap,stap,temp2);
  
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
void hyp_single_rank_dec1_link(su3 link,quad_su3 *conf,int A,int mu,int nu,double alpha1)
{
  if(nrank>0) crash("Works only on single rank!");
  if(mu==nu||nu==rho||mu==rho) crash("The three dirs must be different!");
  
  //find the remaining direction
  int rho=0;
  while(rho==mu || rho==nu) rho++;
  
  //take original link
  su3 temp0;
  su3_prod_double(temp0,conf[A][mu],1-alpha1);
  
  //reset the staple
  su3 stap,temp1,temp2;
  memset(stap,0,sizeof(su3));
  
  //staple in the positive dir
  int B=loclx_neighup[A][eta];
  int F=loclx_neighup[A][mu];
  unsafe_su3_prod_su3(temp1,conf[A][eta],conf[B][mu]);
  unsafe_su3_prod_su3_dag(temp2,temp1,conf[F][eta]);
  su3_summ(stap,stap,temp2);
  
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

void hyp_smear_conf(quad_su3 *smeared_conf,quad_su3 *conf,double alpha0,double alpha1,double alpha2)
{
  communicate_lx_gauge_borders(conf);
  
  //allocate dec1 and dec2
  quad_su3 *dec1_conf,*dec2_conf;
  dec1_conf=nissa_malloc("dec1_conf",loc_vol+loc_bord,quad_su3);
  dec2_conf=nissa_malloc("dec2_conf",loc_vol+loc_bord,quad_su3);
  
  //loop over external index
  for(int mu=0;mu<4;mu++)
    {
      //loop over the first decoration index
      for(int nu=0;nu<4;nu++)
	if(nu!=mu)
	  {
	    //loop over the second decoration index
	    for(int rho=0;rho<4;rho++)
	      if(rho!=mu && rho!=nu)
		{
		  //second level decoration
		  for(int A=0;A<loc_vol;A++)
		    {
		      //take original link
		      su3 temp0;
		      su3_prod_double(temp0,conf[A][mu],1-alpha2);
		      
		      //reset the staple
		      su3 stap,temp1,temp2;
		      memset(stap,0,sizeof(su3));
		      
		      //staple in the positive dir
		      int B=loclx_neighup[A][eta];
		      int F=loclx_neighup[A][mu];
		      unsafe_su3_prod_su3(temp1,conf[A][eta],conf[B][mu]);
		      unsafe_su3_prod_su3_dag(temp2,temp1,conf[F][eta]);
		      su3_summ(stap,stap,temp2);
		      
		      //staple in the negative dir
		      int D=loclx_neighdw[A][eta];
		      int E=loclx_neighup[D][mu];
		      unsafe_su3_dag_prod_su3(temp1,conf[D][eta],conf[D][mu]);
		      unsafe_su3_prod_su3(temp2,temp1,conf[E][eta]);
		      su3_summ(stap,stap,temp2);
		      
		      //summ the two staples with appropriate coef
		      su3_summ_the_prod_double(temp0,stap,alpha2/2);
		      
		      //project the resulting link onto su3
		      su3_unitarize_maximal_trace_projecting(dec2_conf[A][eta],temp0);
		    }
		}
	    
	    //first level decoration
	    for(int A=0;A<loc_vol;A++)
	      {
		//take original link
		su3 temp0;
		su3_prod_double(temp0,conf[A][mu],1-alpha1);
		
		//reset the staple
		su3 stap;
		memset(stap,0,sizeof(su3));
		
		//loop on the rho dir
		for(int rho=0;rho<4;rho++)
		  if(rho!=mu && rho!=nu)
		    {
		      su3 temp1,temp2;
		      
		      //staple in the positive dir
		      int B=loclx_neighup[A][rho];
		      int F=loclx_neighup[A][mu];
		      unsafe_su3_prod_su3(temp1,dec2_conf[A][rho],dec2_conf[B][mu]);
		      unsafe_su3_prod_su3_dag(temp2,temp1,dec2_conf[F][rho]);
		      su3_summ(stap,stap,temp2);
		      
		      //staple in the negative dir
		      int D=loclx_neighdw[A][rho];
		      int E=loclx_neighup[D][mu];
		      unsafe_su3_dag_prod_su3(temp1,dec2_conf[D][rho],dec2_conf[D][mu]);
		      unsafe_su3_prod_su3(temp2,temp1,dec2_conf[E][rho]);
		      su3_summ(stap,stap,temp2);
		      
		      //summ the two staples with appropriate coef
		      su3_summ_the_prod_double(temp0,stap,alpha1/4);
		      
		      //project the resulting link onto su3
		      su3_unitarize_maximal_trace_projecting(dec1_conf[A][rho],temp0);
		    }
	  }
    }
  
  //free dec 1 and dec2
  nissa_free(dec1_conf);
  nissa_free(dec2_conf);
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

  //smeare the conf
  quad_su3 *hyp_smeared_conf=nissa_malloc("hyp_smeared_conf",loc_vol+loc_bord,quad_su3);  
  hyp_smear_conf(hyp_smeared_conf,conf,0.1,0.2,0.3);
  
  nissa_free(conf);
  
  close_nissa();
  
  return 0;
}
