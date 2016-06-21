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

#include "clover_term.hpp"

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
  
  void apply_point_squared_twisted_clover_term_to_halfspincolor(halfspincolor out,double mass,double kappa,clover_term_t Cl,halfspincolor in)
  {
    halfspincolor temp;
    apply_point_twisted_clover_term_to_halfspincolor(temp,+mass,kappa,Cl,in);
    apply_point_twisted_clover_term_to_halfspincolor(out ,-mass,kappa,Cl,temp);
  }
  
  //form the inverse of the clover term
  void invert_point_twisted_clover_term(inv_clover_term_t inv,double mass,double kappa,clover_term_t Cl)
  {
    //inv_clover_term_t dir;
    
    for(int x_high_low=0;x_high_low<2;x_high_low++)
      for(int x_id=0;x_id<NDIRAC/2;x_id++)
	for(int x_ic=0;x_ic<NCOL;x_ic++)
	  {
	    //prepare the point source
	    halfspincolor b;
	    halfspincolor_put_to_zero(b);
	    b[x_id][x_ic][RE]=1;
	    double ori_rr=halfspincolor_norm2(b);
	    
	    //halfspincolor cicc;
	    //apply_point_twisted_clover_term_to_halfspincolor(cicc,mass,kappa,Cl+2*x_high_low,b);
	    //for(int id=0;id<NDIRAC/2;id++)
	    //  for(int ic=0;ic<NCOL;ic++)
	    //	for(int ri=0;ri<2;ri++)
	    //	  dir[x_high_low][id][ic][x_id][x_ic][ri]=cicc[id][ic][ri];
	    
	    //reset the solution
	    halfspincolor x;
	    halfspincolor_put_to_zero(x);
	    
	    //prepare r and p
	    halfspincolor r,p;
	    halfspincolor_copy(r,b);
	    halfspincolor_copy(p,b);
	    
	    //norm of r is 1
	    double rr=1;
	    
	    //count iterations
	    int iter=1;
	    const int niter_max=20;
	    const double target_res=1e-32;
	    do
	      {
		//compute (p,Ap)
		halfspincolor ap;
		apply_point_squared_twisted_clover_term_to_halfspincolor(ap,mass,kappa,Cl+2*x_high_low,p);
		double pap=halfspincolor_scal_prod(p,ap);
		
		//compute alpha, store rr as old one
		double alpha=rr/pap;
		double roro=rr;
		
		//adjust new solution and residual,
		//compute new residual norm
		halfspincolor_summ_the_prod_double(x, p,+alpha);
		halfspincolor_summ_the_prod_double(r,ap,-alpha);
		rr=halfspincolor_norm2(r);
		
		//adjust new krylov vector
		double beta=rr/roro;
		halfspincolor_summ_the_prod_double(p,r,p,beta);
		
		//write residual
		//master_printf("iter %d rel residue: %lg\n",iter,rr/ori_rr);
		iter++;
	      }
	    while(rr/ori_rr>=target_res && iter<niter_max);
	    if(iter>=niter_max) crash("exceeded maximal number of iterations %d, arrived to %d",niter_max,iter);
	    
	    //halfspincolor ap;
	    //apply_point_squared_twisted_clover_term_to_halfspincolor(ap,mass,kappa,Cl+2*x_high_low,x);
	    //halfspincolor_summ_the_prod_double(ap,b,-1);
	    //double trr=halfspincolor_norm2(ap)/ori_rr;
	    //master_printf("true residue: %lg vs %lg\n",trr,rr);
	    
	    //copy the solution after removing the hermitian
	    halfspincolor temp;
	    apply_point_twisted_clover_term_to_halfspincolor(temp,-mass,kappa,Cl+2*x_high_low,x);
	    for(int id=0;id<NDIRAC/2;id++)
	      for(int ic=0;ic<NCOL;ic++)
		for(int ri=0;ri<2;ri++)
		  inv[x_high_low][id][ic][x_id][x_ic][ri]=temp[id][ic][ri];
	  }
    
    // for(int x_high_low=0;x_high_low<2;x_high_low++)
    //   {
    // 	for(int x_id=0;x_id<NDIRAC/2;x_id++)
    // 	  for(int x_ic=0;x_ic<NCOL;x_ic++)
    // 	    {
    // 	      for(int id=0;id<NDIRAC/2;id++)
    // 		for(int ic=0;ic<NCOL;ic++)
    // 		  {
    // 		    for(int ri=0;ri<2;ri++)
    // 		      master_printf("%+2.2g ",inv[x_high_low][x_id][x_ic][id][ic][ri]);
    // 		    master_printf(" ");
    // 		  }
    // 	      master_printf("\n\n");
    // 	    }
    // 	master_printf("\n\n");
    //   }
    
    //for(int x_high_low=0;x_high_low<2;x_high_low++)
    // {
    //	for(int x_id=0;x_id<NDIRAC/2;x_id++)
    //	  for(int x_ic=0;x_ic<NCOL;x_ic++)
    //	    {
    //	      for(int id=0;id<NDIRAC/2;id++)
    //		for(int ic=0;ic<NCOL;ic++)
    //		  {
    //		    for(int ri=0;ri<2;ri++)
    //		      master_printf("%+2.2g ",dir[x_high_low][x_id][x_ic][id][ic][ri]);
    //		    master_printf(" ");
    //		  }
    //	      master_printf("\n\n");
    //	    }
    //	master_printf("\n\n");
    //  }
  }
  
  THREADABLE_FUNCTION_4ARG(invert_twisted_clover_term, inv_clover_term_t*,invCl, double,mass, double,kappa, clover_term_t*,Cl)
  {
    GET_THREAD_ID();
    if(IS_MASTER_THREAD) verbosity_lv2_master_printf("Computing inverse clover term for quark of mass %lg and kappa %lg\n",mass,kappa);
    NISSA_PARALLEL_LOOP(X,0,loc_vol) invert_point_twisted_clover_term(invCl[X],mass,kappa,Cl[X]);
    set_borders_invalid(invCl);
  }
  THREADABLE_FUNCTION_END

}
