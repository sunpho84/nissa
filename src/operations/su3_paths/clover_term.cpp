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
  CUDA_HOST_AND_DEVICE void point_chromo_operator(clover_term_t Cl,quad_su3 *conf,int X)
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
  void chromo_operator(clover_term_t* Cl,quad_su3* conf)
  {
    master_printf("Computing Chromo operator\n");
    communicate_lx_quad_su3_edges(conf);
    NISSA_PARALLEL_LOOP(X,0,locVol)
      point_chromo_operator(Cl[X],conf,X);
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(Cl);
  }
  
  void chromo_operator(eo_ptr<clover_term_t> Cl_eo,eo_ptr<quad_su3> conf_eo)
  {
    quad_su3 *conf_lx=nissa_malloc("conf_lx",locVol+bord_vol+edge_vol,quad_su3);
    clover_term_t *Cl_lx=nissa_malloc("Cl_lx",locVol,clover_term_t);
    paste_eo_parts_into_lx_vector(conf_lx,conf_eo);
    chromo_operator(Cl_lx,conf_lx);
    split_lx_vector_into_eo_parts(Cl_eo,Cl_lx);
    nissa_free(conf_lx);
    nissa_free(Cl_lx);
  }
  
  //apply the chromo operator to the passed spincolor
  CUDA_HOST_AND_DEVICE void unsafe_apply_point_chromo_operator_to_spincolor(spincolor out,clover_term_t Cl,spincolor in)
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
  void unsafe_apply_chromo_operator_to_spincolor(spincolor* out,clover_term_t* Cl,spincolor* in)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      unsafe_apply_point_chromo_operator_to_spincolor(out[ivol],Cl[ivol],in[ivol]);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
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
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      unsafe_apply_point_chromo_operator_to_spincolor_128(out[ivol],Cl[ivol],in[ivol]);
    NISSA_PARALLEL_LOOP_END;
    set_borders_invalid(out);
  }
  
  //apply the chromo operator to the passed colorspinspin
  //normalization as in ape next
  void unsafe_apply_chromo_operator_to_colorspinspin(colorspinspin* out,clover_term_t* Cl,colorspinspin* in)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	spincolor temp1,temp2;
	
	//Loop over the four source dirac indexes
	for(int id_source=0;id_source<NDIRAC;id_source++) //dirac index of source
	  {
	    //Switch the color_spinspin into the spincolor.
	    get_spincolor_from_colorspinspin(temp1,in[ivol],id_source);
	    
	    unsafe_apply_point_chromo_operator_to_spincolor(temp2,Cl[ivol],temp1);
	    
	    //Switch back the spincolor into the colorspinspin
	    put_spincolor_into_colorspinspin(out[ivol],temp2,id_source);
	  }
      }
    NISSA_PARALLEL_LOOP_END;
    
    //invalidate borders
    set_borders_invalid(out);
  }
  
  //apply the chromo operator to the passed su3spinspin
  //normalization as in ape next
  void unsafe_apply_chromo_operator_to_su3spinspin(su3spinspin* out,clover_term_t* Cl,su3spinspin* in)
  {
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	spincolor temp1,temp2;
	
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
      }
    NISSA_PARALLEL_LOOP_END;
    
    //invalidate borders
    set_borders_invalid(out);
  }
  
  //apply a diagonal matrix plus clover term to up or low components
  CUDA_HOST_AND_DEVICE void apply_point_diag_plus_clover_term_to_halfspincolor(halfspincolor out,complex diag,clover_term_t Cl,halfspincolor in)
  {
    unsafe_color_prod_complex(out[0],in[0],diag);
    su3_summ_the_prod_color(out[0],Cl[0],in[0]);
    su3_dag_summ_the_prod_color(out[0],Cl[1],in[1]);
    
    unsafe_color_prod_complex(out[1],in[1],diag);
    su3_summ_the_prod_color(out[1],Cl[1],in[0]);
    su3_subt_the_prod_color(out[1],Cl[0],in[1]);
  }
  CUDA_HOST_AND_DEVICE void apply_point_diag_plus_clover_term_to_halfspincolor_128(halfspincolor_128 out,complex& diag,clover_term_t Cl,halfspincolor_128 in)
  {
    unsafe_color_128_prod_complex_64(out[0],in[0],diag);
    su3_summ_the_prod_color_128(out[0],Cl[0],in[0]);
    su3_dag_summ_the_prod_color_128(out[0],Cl[1],in[1]);
    
    unsafe_color_128_prod_complex_64(out[1],in[1],diag);
    su3_summ_the_prod_color_128(out[1],Cl[1],in[0]);
    su3_subt_the_prod_color_128(out[1],Cl[0],in[1]);
  }
  
  CUDA_HOST_AND_DEVICE void apply_point_squared_twisted_clover_term_to_halfspincolor(halfspincolor out,double mass,double kappa,clover_term_t Cl,halfspincolor in)
  {
    halfspincolor temp;
    apply_point_twisted_clover_term_to_halfspincolor(temp,+mass,kappa,Cl,in);
    apply_point_twisted_clover_term_to_halfspincolor(out ,-mass,kappa,Cl,temp);
  }
  
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
    // master_printf("DIFF: %lg\n",sqrt(diff/norm2));
    
    // auto pr=[](halfspincolor_halfspincolor out)
    // 	    {
    // 	      for(int id1=0;id1<NDIRAC/2;id1++)
    // 		for(int ic1=0;ic1<NCOL;ic1++)
    // 		  {
    // 		    for(int id=0;id<NDIRAC/2;id++)
    // 		      for(int ic=0;ic<NCOL;ic++)
    // 			master_printf("%lg %lg\t",out[id1][ic1][id][ic][RE],out[id1][ic1][id][ic][IM]);
    // 		    master_printf("\n");
    // 		  }
    // 		    master_printf("\n");
    // 	    };
    
    // master_printf("out:\n");
    // pr(out);
    // master_printf("correct:\n");
    // pr(out_sure);
  }
  
  //form the inverse of the clover term
  CUDA_HOST_AND_DEVICE void invert_point_twisted_clover_term(inv_clover_term_t inv,double mass,double kappa,clover_term_t Cl)
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
	    [[ maybe_unused ]]
	    const int niter_max=200,niter_for_verbosity=20;
	    const double target_res=1e-32;
	    double res;
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
		
		//compute resiude
		res=rr/ori_rr;
		
		//write residual
		if(iter>=niter_for_verbosity)
#ifndef COMPILING_FOR_DEVICE
		  master_printf("iter %d rel residue: %lg\n",iter,res)
#endif
		    ;
		iter++;
	      }
	    while(res>=target_res && iter<niter_max);
	    if(iter>=niter_max)
#ifndef COMPILING_FOR_DEVICE
	      crash("exceeded maximal number of iterations %d, arrived to %d with residue %lg, target %lg",niter_max,iter,res,target_res);
 #else
	    __trap();
#endif
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
}

// #include "dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec.hpp"
// #include "geometry/geometry_mix.hpp"

namespace nissa
 {
   void invert_twisted_clover_term(inv_clover_term_t* invCl,double mass,double kappa,clover_term_t* Cl)
   {
     if(IS_MASTER_THREAD) verbosity_lv2_master_printf("Computing inverse clover term for quark of mass %lg and kappa %lg\n",mass,kappa);
     NISSA_PARALLEL_LOOP(X,0,get_vect(invCl)->nel)
       invert_point_twisted_clover_term(invCl[X],mass,kappa,Cl[X]);
     NISSA_PARALLEL_LOOP_END;
     
     set_borders_invalid(invCl);
     
     // //check
     // inv_clover_term_t *invCl_eo[2]={nissa_malloc("inver",loc_volh+bord_volh,inv_clover_term_t),nissa_malloc("inver",loc_volh+bord_volh,inv_clover_term_t)};
     // clover_term_t *Cl_eo[2]={nissa_malloc("Cleo",loc_volh+bord_volh,clover_term_t),nissa_malloc("Cleo",loc_volh+bord_volh,clover_term_t)};
     // spincolor *source[2]={nissa_malloc("source",loc_volh+bord_volh,spincolor),nissa_malloc("source",loc_volh+bord_volh,spincolor)};
     // generate_fully_undiluted_eo_source(source,RND_Z4,-1);
     // spincolor *inver[2]={nissa_malloc("inver",loc_volh+bord_volh,spincolor),nissa_malloc("inver",loc_volh+bord_volh,spincolor)};
     // split_lx_vector_into_eo_parts(invCl_eo,invCl);
     // split_lx_vector_into_eo_parts(Cl_eo,Cl);
     // spincolor *reco[2]={nissa_malloc("reco",loc_volh+bord_volh,spincolor),nissa_malloc("reco",loc_volh+bord_volh,spincolor)};
     // for(int eo=0;eo<2;eo++)
     //   {
     // 	 inv_tmclovDee_or_oo_eos(inver[eo],invCl_eo[eo],false,source[eo]);
     // 	 tmclovDee_or_oo_eos(reco[eo],kappa,Cl_eo[eo],false,mass,inver[eo]);
     // 	 double_vector_subtassign((double*)(reco[eo]),(double*)(source[eo]),loc_volh*sizeof(spincolor)/sizeof(double));
     // 	 double r=double_vector_glb_norm2(reco[eo],loc_volh);
     // 	 master_printf("Check of the inverse of Clover: %lg\n",r);
     //   }
     
     // for(int eo=0;eo<2;eo++)
     //   {
     // 	 nissa_free(source[eo]);
     // 	 nissa_free(inver[eo]);
     // 	 nissa_free(reco[eo]);
     // 	 nissa_free(Cl_eo[eo]);
     // 	 nissa_free(invCl_eo[eo]);
     //   }
   }
}
