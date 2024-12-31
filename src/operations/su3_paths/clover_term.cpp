#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include "base/vectors.hpp"
#include "communicate/edges.hpp"
#include "linalgs/linalgs.hpp"
#include "measures/gauge/topological_charge.hpp"
#include "new_types/float_128.hpp"

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
  template <typename T,
	    typename U>
  CUDA_HOST_AND_DEVICE
  void point_chromo_operator(T&& Cl,
			     const U& conf,
			     const int& X)
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
	    P[i][ic1][ic2][RE]=(leaves[i][ic1][ic2][RE]-leaves[i][ic2][ic1][RE])/4;
	    P[i][ic1][ic2][IM]=(leaves[i][ic1][ic2][IM]+leaves[i][ic2][ic1][IM])/4;
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
  
  void chromo_operator(LxField<clover_term_t>& Cl,
		       const LxField<quad_su3>& conf)
  {
    MASTER_PRINTF("Computing Chromo operator\n");
    
    conf.updateEdges();
    PAR(0,locVol,
	CAPTURE(TO_READ(conf),
		TO_WRITE(Cl)),
	X,
	{
	  point_chromo_operator(Cl[X],conf,X);
	});
  }
  
  void chromo_operator(EoField<clover_term_t>& Cl_eo,
		       const EoField<quad_su3>& conf_eo)
  {
    LxField<quad_su3> conf_lx("conf_lx",WITH_HALO_EDGES);
    LxField<clover_term_t> Cl_lx("Cl_lx");
    
    paste_eo_parts_into_lx_vector(conf_lx,conf_eo);
    chromo_operator(Cl_lx,conf_lx);
    
    split_lx_vector_into_eo_parts(Cl_eo,Cl_lx);
  }
  
  void unsafe_apply_chromo_operator_to_spincolor(LxField<spincolor>& out,
						 const LxField<clover_term_t>& Cl,
						 const LxField<spincolor>& in)
  {
    PAR(0,locVol,
	CAPTURE(TO_WRITE(out),
		TO_READ(Cl),
		TO_READ(in)),
	ivol,
	{
	  unsafe_apply_point_chromo_operator_to_spincolor(out[ivol],Cl[ivol],in[ivol]);
	});
  }
  
  //128 bit case
  CUDA_HOST_AND_DEVICE void unsafe_apply_point_chromo_operator_to_spincolor_128(spincolor_128 out,clover_term_t Cl,spincolor_128 in)
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
  void unsafe_apply_chromo_operator_to_spincolor_128(spincolor_128* out,clover_term_t* Cl,spincolor_128* in)
  {
    CRASH("reimplement");
    
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   unsafe_apply_point_chromo_operator_to_spincolor_128(out[ivol],Cl[ivol],in[ivol]);
    // NISSA_PARALLEL_LOOP_END;
    // set_borders_invalid(out);
  }
  
  //apply the chromo operator to the passed colorspinspin
  //normalization as in ape next
  void unsafe_apply_chromo_operator_to_colorspinspin(colorspinspin* out,clover_term_t* Cl,colorspinspin* in)
  {
    CRASH("reimplement");
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   {
    // 	spincolor temp1,temp2;
	
    // 	//Loop over the four source dirac indexes
    // 	for(int id_source=0;id_source<NDIRAC;id_source++) //dirac index of source
    // 	  {
    // 	    //Switch the color_spinspin into the spincolor.
    // 	    get_spincolor_from_colorspinspin(temp1,in[ivol],id_source);
	    
    // 	    unsafe_apply_point_chromo_operator_to_spincolor(temp2,Cl[ivol],temp1);
	    
    // 	    //Switch back the spincolor into the colorspinspin
    // 	    put_spincolor_into_colorspinspin(out[ivol],temp2,id_source);
    // 	  }
    //   }
    // NISSA_PARALLEL_LOOP_END;
    
    // //invalidate borders
    // set_borders_invalid(out);
  }
  
  //apply the chromo operator to the passed su3spinspin
  //normalization as in ape next
  void unsafe_apply_chromo_operator_to_su3spinspin(su3spinspin* out,clover_term_t* Cl,su3spinspin* in)
  {
    CRASH("reimplement");
    
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   {
    // 	spincolor temp1,temp2;
	
    // 	//Loop over the four source dirac indexes
    // 	for(int id_source=0;id_source<NDIRAC;id_source++) //dirac index of source
    // 	  for(int ic_source=0;ic_source<NCOL;ic_source++) //color index of source
    // 	    {
    // 	      //Switch the su3spinspin into the spincolor.
    // 	      get_spincolor_from_su3spinspin(temp1,in[ivol],id_source,ic_source);
	      
    // 	      unsafe_apply_point_chromo_operator_to_spincolor(temp2,Cl[ivol],temp1);
	      
    // 	      //Switch back the spincolor into the colorspinspin
    // 	      put_spincolor_into_su3spinspin(out[ivol],temp2,id_source,ic_source);
    // 	    }
    //   }
    // NISSA_PARALLEL_LOOP_END;
    
    // //invalidate borders
    // set_borders_invalid(out);
  }
  
  // CUDA_HOST_AND_DEVICE void apply_point_diag_plus_clover_term_to_halfspincolor_128(halfspincolor_128 out,complex& diag,clover_term_t Cl,halfspincolor_128 in)
  // {
  //   unsafe_color_128_prod_complex_64(out[0],in[0],diag);
  //   su3_summ_the_prod_color_128(out[0],Cl[0],in[0]);
  //   su3_dag_summ_the_prod_color_128(out[0],Cl[1],in[1]);
    
  //   unsafe_color_128_prod_complex_64(out[1],in[1],diag);
  //   su3_summ_the_prod_color_128(out[1],Cl[1],in[0]);
  //   su3_subt_the_prod_color_128(out[1],Cl[0],in[1]);
  // }
  
  CUDA_HOST_AND_DEVICE void fill_point_twisted_clover_term(halfspincolor_halfspincolor out,int x_high_low,clover_term_t C,double mass,double kappa)
  {
    // halfspincolor_halfspincolor out_sure;
    // for(int id1=0;id1<NDIRAC/2;id1++)
    //   for(int ic1=0;ic1<NCOL;ic1++)
    // 	{
    // 	  halfspincolor in;
    // 	  halfspincolor_put_to_zero(in);
    // 	  in[id1][ic1][RE]=1.0;
	  
    // 	  halfspincolor temp;
    // 	  apply_point_twisted_clover_term_to_halfspincolor(temp,mass,kappa,C+2*x_high_low,in);
	  
    // 	  for(int id=0;id<NDIRAC/2;id++)
    // 	    for(int ic=0;ic<NCOL;ic++)
    // 	      complex_copy(out_sure[id][ic][id1][ic1],temp[id][ic]);
    // 	}
    
    for(int ic1=0;ic1<NCOL;ic1++)
      for(int ic2=0;ic2<NCOL;ic2++)
	{
	  auto fill=[&out,x_high_low,icl1=ic1,icl2=ic2,C](int id1,int id2,int icl,int ic1,int ic2,const complex& f)
		    {
		      for(int ri=0;ri<2;ri++)
			out[id1][icl1][id2][icl2][ri]=C[icl+x_high_low*2][ic1][ic2][ri]*f[ri];
		    };
	  
	  fill(0,0, 0,ic1,ic2, {1,1});
	  fill(1,0, 1,ic1,ic2, {1,1});
	  fill(0,1, 1,ic2,ic1, {1,-1});
	  fill(1,1, 0,ic1,ic2, {-1,-1});
	}
    
    complex mt={1/(2*kappa),(x_high_low==0)?mass:-mass};
    for(int id=0;id<NDIRAC/2;id++)
      for(int ic=0;ic<NCOL;ic++)
	complex_summassign(out[id][ic][id][ic],mt);
    
    // double diff=0,norm2=0;
    // for(int id1=0;id1<NDIRAC/2;id1++)
    //   for(int ic1=0;ic1<NCOL;ic1++)
    // 	for(int id=0;id<NDIRAC/2;id++)
    // 	  for(int ic=0;ic<NCOL;ic++)
    // 	    {
    // 	      norm2+=complex_norm2(out_sure[id][ic][id1][ic1]);
    // 	      complex_subtassign(out_sure[id][ic][id1][ic1],out[id][ic][id1][ic1]);
    // 	      diff+=complex_norm2(out_sure[id][ic][id1][ic1]);
    // 	    }
    // MASTER_PRINTF("DIFF: %lg\n",sqrt(diff/norm2));
    
    // auto pr=[](halfspincolor_halfspincolor out)
    // 	    {
    // 	      for(int id1=0;id1<NDIRAC/2;id1++)
    // 		for(int ic1=0;ic1<NCOL;ic1++)
    // 		  {
    // 		    for(int id=0;id<NDIRAC/2;id++)
    // 		      for(int ic=0;ic<NCOL;ic++)
    // 			MASTER_PRINTF("%lg %lg\t",out[id1][ic1][id][ic][RE],out[id1][ic1][id][ic][IM]);
    // 		    MASTER_PRINTF("\n");
    // 		  }
    // 		    MASTER_PRINTF("\n");
    // 	    };
    
    // MASTER_PRINTF("out:\n");
    // pr(out);
    // MASTER_PRINTF("correct:\n");
    // pr(out_sure);
  }
}
