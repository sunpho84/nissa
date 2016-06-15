#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/edges.hpp"
#include "linalgs/linalgs.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "topological_charge.hpp"

namespace nissa
{
  /*
    
    build the chromo operator
    
    i \sigma_\mu_\nu F_\mu_\nu
    
    This is organised as follow. The six independent component of
    P_\mu_\nu are stored in P_i, with each "i" associated to the 6 independent
    gammas (10-15),
    
    from anti-symmetrized four-leaves, following this decomposition
    
    
    (0,-P2+P3),     (-P1-P4,-P0+P5), (0,0),          (0,0)
    (P1+P4,-P0+P5), (0,P2-P3),       (0,0),          (0,0)
    (0,0),          (0,0),           (0,P2+P3),      (P1-P4,P0+P5)
    (0,0),          (0,0),           (-P1+P4,P0+P5), (0,-P2-P3)
    
    +iA   -B+iC  0      0
    +B+iC -iA    0      0
    0      0     +iD    -E+iF
    0      0     +E+iF  -iD
    
    A=P3-P2, B=P4+P1, C=P5-P0
    D=P3+P2, E=P4-P1, F=P5+P0
    
    +G  +H^+ 0    0
    +H  -G   0    0
    0    0   +I  +J^+
    0    0   +J  -I
    
    out[0]=G=iA, out[1]=H=B+iC
    out[2]=I=iD, out[3]=J=E+iF
    
    NB: indeed Pi is anti-hermitian
  */
  void point_chromo_operator(clover_term_t Cl,quad_su3 *conf,int X)
  {
    //this is the non-anti-symmetric part 2*F_mu_nu
    as2t_su3 leaves;
    four_leaves_point(leaves,conf,X);
    
    //anti-symmetrize and divide by 2
    as2t_su3 P;
    for(int i=0;i<6;i++)
      for(int ic1=0;ic1<NCOL;ic1++)
	for(int ic2=0;ic2<NCOL;ic2++)
	  {
	    P[i][ic1][ic2][0]=(leaves[i][ic1][ic2][0]-leaves[i][ic2][ic1][0])/4;
	    P[i][ic1][ic2][1]=(leaves[i][ic1][ic2][1]+leaves[i][ic2][ic1][1])/4;
	  }
    
    su3 A,B,C,D,E,F;
    su3_subt(A,P[3],P[2]);
    su3_summ(B,P[4],P[1]);
    su3_subt(C,P[5],P[0]);
    su3_summ(D,P[3],P[2]);
    su3_subt(E,P[4],P[1]);
    su3_summ(F,P[5],P[0]);
    
    //0
    su3_prod_idouble(Cl[0],A,1);
    //1
    su3_copy(Cl[1],B);
    su3_summ_the_prod_idouble(Cl[1],C,1);
    //2
    su3_prod_idouble(Cl[2],D,1);
    //3
    su3_copy(Cl[3],E);
    su3_summ_the_prod_idouble(Cl[3],F,1);
    
  }
  THREADABLE_FUNCTION_2ARG(chromo_operator, clover_term_t*,Cl, quad_su3*,conf)
  {
    GET_THREAD_ID();
    communicate_lx_quad_su3_edges(conf);
    NISSA_PARALLEL_LOOP(X,0,loc_vol) point_chromo_operator(Cl[X],conf,X);
    set_borders_invalid(Cl);
  }
  THREADABLE_FUNCTION_END
  
  //apply the chromo operator to the passed spincolor
  void unsafe_apply_point_chromo_operator_to_spincolor(spincolor out,clover_term_t Cl,spincolor in)
  {
    unsafe_su3_prod_color(out[0],Cl[0],in[0]);
    su3_dag_summ_the_prod_color(out[0],Cl[1],in[1]);
    unsafe_su3_prod_color(out[1],Cl[1],in[0]);
    su3_subt_the_prod_color(out[1],Cl[0],in[1]);
    
    unsafe_su3_prod_color(out[2],Cl[2],in[2]);
    su3_dag_summ_the_prod_color(out[2],Cl[3],in[3]);
    unsafe_su3_prod_color(out[3],Cl[3],in[2]);
    su3_subt_the_prod_color(out[3],Cl[2],in[3]);
  }
  THREADABLE_FUNCTION_3ARG(unsafe_apply_chromo_operator_to_spincolor, spincolor*,out, clover_term_t*,Cl, spincolor*,in)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      unsafe_apply_point_chromo_operator_to_spincolor(out[ivol],Cl[ivol],in[ivol]);
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //128 bit case
  void unsafe_apply_point_chromo_operator_to_spincolor_128(spincolor_128 out,clover_term_t Cl,spincolor_128 in)
  {
    unsafe_su3_prod_color_128(out[0],Cl[0],in[0]);
    su3_dag_summ_the_prod_color_128(out[0],Cl[1],in[1]);
    unsafe_su3_prod_color_128(out[1],Cl[1],in[0]);
    su3_subt_the_prod_color_128(out[1],Cl[0],in[1]);
    
    unsafe_su3_prod_color_128(out[2],Cl[2],in[2]);
    su3_dag_summ_the_prod_color_128(out[2],Cl[3],in[3]);
    unsafe_su3_prod_color_128(out[3],Cl[3],in[2]);
    su3_subt_the_prod_color_128(out[3],Cl[2],in[3]);
  }
  THREADABLE_FUNCTION_3ARG(unsafe_apply_chromo_operator_to_spincolor_128, spincolor_128*,out, clover_term_t*,Cl, spincolor_128*,in)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      unsafe_apply_point_chromo_operator_to_spincolor_128(out[ivol],Cl[ivol],in[ivol]);
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //apply the chromo operator to the passed colorspinspin
  //normalization as in ape next
  THREADABLE_FUNCTION_3ARG(unsafe_apply_chromo_operator_to_colorspinspin, colorspinspin*,out, clover_term_t*,Cl, colorspinspin*,in)
  {
    spincolor temp1,temp2;
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      //Loop over the four source dirac indexes
      for(int id_source=0;id_source<NDIRAC;id_source++) //dirac index of source
	{
	  //Switch the color_spinspin into the spincolor.
	  get_spincolor_from_colorspinspin(temp1,in[ivol],id_source);
	  
	  unsafe_apply_point_chromo_operator_to_spincolor(temp2,Cl[ivol],temp1);
	  
	  //Switch back the spincolor into the colorspinspin
	  put_spincolor_into_colorspinspin(out[ivol],temp2,id_source);
	}
    
    //invalidate borders
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //apply the chromo operator to the passed su3spinspin
  //normalization as in ape next
  THREADABLE_FUNCTION_3ARG(unsafe_apply_chromo_operator_to_su3spinspin, su3spinspin*,out, clover_term_t*,Cl, su3spinspin*,in)
  {
    spincolor temp1,temp2;
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      //Loop over the four source dirac indexes
      for(int id_source=0;id_source<NDIRAC;id_source++) //dirac index of source
	for(int ic_source=0;ic_source<NCOL;ic_source++) //color index of source
	  {
	    //Switch the su3spinspin into the spincolor.
	    get_spincolor_from_su3spinspin(temp1,in[ivol],id_source,ic_source);
	    
	    unsafe_apply_point_chromo_operator_to_spincolor(temp2,Cl[ivol],temp1);
	    
	    //Switch back the spincolor into the colorspinspin
	    put_spincolor_into_su3spinspin(out[ivol],temp2,id_source,ic_source);
	  }
    
    //invalidate borders
    set_borders_invalid(out);
  }
  THREADABLE_FUNCTION_END
  
  //apply a diagonal matrix plus clover term to up or low components
  void apply_point_diag_plus_clover_term_to_halfspincolor(halfspincolor out,complex diag,clover_term_t Cl,halfspincolor in)
  {
    unsafe_color_prod_complex(out[0],in[0],diag);
    su3_summ_the_prod_color(out[0],Cl[0],in[0]);
    su3_dag_summ_the_prod_color(out[0],Cl[1],in[1]);
    
    unsafe_color_prod_complex(out[1],in[1],diag);
    su3_summ_the_prod_color(out[1],Cl[1],in[0]);
    su3_subt_the_prod_color(out[1],Cl[0],in[1]);
  }
  
  void apply_point_twisted_clover_term_to_halfspincolor(halfspincolor out,double mass,double kappa,clover_term_t Cl,halfspincolor in)
  {
    
  }
  
  void apply_point_squared_twisted_clover_term_to_halfspincolor(halfspincolor out,double mass,double kappa,clover_term_t Cl,halfspincolor in)
  {
    
  }
  
  double halfspincolor_scal_prod(halfspincolor a,halfspincolor b)
  {
    double out=0;
    for(int id=0;id<NDIRAC/2;id++)
      for(int icol=0;icol<NCOL;icol++)
	for(int ri=0;ri<2;ri++)
	  out+=a[id][icol][ri]*b[id][icol][ri];
	    
    return out;
  }
  
  void halfspincolor_copy(halfspincolor a,halfspincolor b)
  {memcpy(a,b,sizeof(halfspincolor));}
  
  //form the inverse of the clover term
  void invert_point_twisted_clover_term(oct_su3 inv,double mass,double kappa,double cSW,clover_term_t Cl)
  {
    int nd=NDIRAC/2*NCOL*2;
    
    for(int high_low=0;high_low<2;high_low++)
      for(int x_id=0;x_id<NDIRAC/2;x_id++)
	for(int x_ic=0;x_ic<NCOL;x_ic++)
	  {
	    //prepare the source
	    double b[nd];
	    memset(b,0,nd*sizeof(double));
	    b[RE+2*(x_ic+NCOL*x_id)]=1;
	    
	    //reset the solution
	    double x[nd];
	    memset(x,0,nd*sizeof(double));
	    
	    //prepare r and p
	    double r[nd],p[nd];
	    memcpy(r,b,nd*sizeof(double));
	    memcpy(p,b,nd*sizeof(double));
	    
	    //norm of r is 1
	    double rr=1;
	    
	    //count iterations
	    int iter=1;
	    do
	      {
		//compute (p,Ap)
		double ap[nd];
		apply_point_twisted_clover_term_to_halfspincolor((color*)ap,mass,kappa,Cl,(color*)p);
		double pap=0;
		for(int i=0;i<nd;i++) pap+=p[i]*ap[i];
		
		//compute alpha, store rr as old one
		double alpha=rr/pap;
		double roro=rr;
		
		//adjust new solution and residual,
		//compute new residual norm
		rr=0;
		for(int i=0;i<nd;i++)
		  {
		    x[i]+=alpha*p[i];
		    r[i]-=alpha*ap[i];
		    rr+=r[i]*r[i]; //computing residual here we save one loop
		  }
		
		//adjust new krylov vector
		double beta=rr/roro;
		for(int i=0;i<nd;i++) p[i]=r[i]+beta*p[i];
		
		//write residual
		master_printf("iter %d residue: %lg\n",iter,rr);
		iter++;
	      }
	    while(rr>=1e-32);
	    
	    //copy the solution
	    //out[low]
	  }
  }
  
}
