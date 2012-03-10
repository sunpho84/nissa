#include "nissa.h"

int main(int narg,char **arg)
{
  init_nissa();
  
  //set the lattice grid
  int T=8;
  int L=4;
  
  //init the grid
  init_grid(T,L);
  
  quad_su3 *conf=nissa_malloc("conf",loc_vol+loc_bord+loc_edge,quad_su3);
  
  //read conf, compute plaquette, print it
  read_ildg_gauge_conf(conf,"../../data/L4T8conf");

  int loc_site=0;
  int mu=1;
  double alpha=0.1;

  //calculate staples
  su3 stap,temp1,temp2;
  memset(stap,0,sizeof(su3));
  for(int nu=1;nu<4;nu++)                   //  E---F---C   
    if(nu!=mu)                              //  |   |   | mu
      {                                     //  D---A---B   
	int A=loc_site;                     //        nu    
	int B=loclx_neighup[A][nu];
	int F=loclx_neighup[A][mu];
	unsafe_su3_prod_su3(temp1,conf[A][nu],conf[B][mu]);
	unsafe_su3_prod_su3_dag(temp2,temp1,conf[F][nu]);
	su3_summ(stap,stap,temp2);
        
	int D=loclx_neighdw[A][nu];
	int E=loclx_neighup[D][mu];
	unsafe_su3_dag_prod_su3(temp1,conf[D][nu],conf[D][mu]);
	unsafe_su3_prod_su3(temp2,temp1,conf[E][nu]);
	su3_summ(stap,stap,temp2);
      }
  
  //create new link to be reunitarized
  su3 prop_link;
  for(int icol1=0;icol1<3;icol1++)
    for(int icol2=0;icol2<3;icol2++)
      for(int ri=0;ri<2;ri++)
	prop_link[icol1][icol2][ri]=conf[loc_site][mu][icol1][icol2][ri]+alpha*stap[icol1][icol2][ri];
  
  su3 g;
  su3_unitarize_maximal_trace_projecting(g,prop_link);
  
  nissa_free(conf);
  
  close_nissa();
  
  return 0;
}
