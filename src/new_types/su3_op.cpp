#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "su3_op.hpp"
#include "geometry/geometry_eo.hpp"
#include "operations/su3_paths/rectangular_staples.hpp"

namespace nissa
{
  //gell-mann matrices as from eq.A.10 of Gattringer - note that T=lambda/2
  su3 gell_mann_matr[NCOL*NCOL-1]={
    {
	{{0,0},{1,0},{0,0}},
	{{1,0},{0,0},{0,0}},
	{{0,0},{0,0},{0,0}}},
     {
       {{0,0},{0,-1},{0,0}},
       {{0,1},{0,0},{0,0}},
       {{0,0},{0,0},{0,0}}},
     {
       {{1,0},{0,0},{0,0}},
       {{0,0},{-1,0},{0,0}},
       {{0,0},{0,0},{0,0}}},
     {
       {{0,0},{0,0},{1,0}},
       {{0,0},{0,0},{0,0}},
       {{1,0},{0,0},{0,0}}},
    {
	{{0,0},{0,0},{0,-1}},
	{{0,0},{0,0},{0,0}},
	{{0,1},{0,0},{0,0}}},
     {
       {{0,0},{0,0},{0,0}},
       {{0,0},{0,0},{1,0}},
       {{0,0},{1,0},{0,0}}},
     {
       {{0,0},{0,0},{0,0}},
       {{0,0},{0,0},{0,-1}},
       {{0,0},{0,1},{0,0}}},
     {
       {{1/sqrt(3),0},{0,0},{0,0}},
       {{0,0},{1/sqrt(3),0},{0,0}},
       {{0,0},{0,0},{-2/sqrt(3),0}}}};
  
  CUDA_MANAGED int su3_sub_gr_indices[3][2]={{0,1},{1,2},{0,2}};
  
  //make unitary maximazing Trace(out*M^dag)
  CUDA_HOST_AND_DEVICE void su3_unitarize_maximal_trace_projecting(su3 out,const su3 M,const double precision,const int niter_max)
  {
    //initialize the guess with the identity - proved to be faster than any good guess,
    //because iterations are so good
    su3 U;
    su3_put_to_id(U);
    
    //compute the product
    su3 prod;
    unsafe_su3_prod_su3_dag(prod,U,M);
    
    int iter=0;
    double rotating_norm=1e300;
    int converged;
    do
      {
	converged=true;
	
	for(int overrelax=0;overrelax<3;overrelax++)
	  for(int isub_gr=0;isub_gr<NCOL;isub_gr++)
	    {
	      //take the subgroup isub_gr
	      double r0,r1,r2,r3;
	      su2_part_of_su3(r0,r1,r2,r3,prod,isub_gr);
	      
	      //form the matrix
	      double x0,x1,x2,x3;
	      if(overrelax) su2_get_overrelaxing(x0,x1,x2,x3, r0,r1,r2,r3);
	      else          su2_inv(x0,x1,x2,x3, r0,r1,r2,r3);
	      
	      //modify the subgroup and the product
	      su2_prodassign_su3(x0,x1,x2,x3,isub_gr,U);
	      su2_prodassign_su3(x0,x1,x2,x3,isub_gr,prod);
	      
	      //condition to exit
	      if(!overrelax) rotating_norm=sqrt(su2_nonunitarity(x0,x1,x2,x3));
	      converged&=(iter>=3 and rotating_norm<precision);
	    }
	iter++;
	
	//refix halfway to the end
	double non_un_tol=1e-15;
	double non_un=su3_get_non_unitariness(U);
	if(rotating_norm<sqrt(precision) and non_un>non_un_tol)
	  {
	    su3_unitarize_explicitly_inverting(U,U);
	    unsafe_su3_prod_su3_dag(prod,U,M);
	  }
	
  if(iter>niter_max*0.9)
	  {
	    printf("We arrived to %d iter, that was set to be the maximum\n",iter);
	    printf("Here you are the input link:\n");
	    //su3_print(M);
	    printf("Here you are the current maxtrace link:\n");
	    //su3_print(U);
	    printf("This is meant to be the product:\n");
	    //su3_print(prod);
	    printf("The norm was: %16.16lg and the trace: %16.16lg\n",rotating_norm,su3_real_trace(prod));
	    if(iter>niter_max)
#ifdef COMPILING_FOR_DEVICE
	      __trap();
#else
	      CRASH("%lg",rotating_norm);
#endif
	  }
      }
    while(not converged);
    
    su3_copy(out,U);
  }
  
  //unitarize returning (VV^\dagger)^(-1/2)*V hep-lat/0610092
  void su3_unitarize_with_sqrt(su3 out,const su3 in)
  {
#ifdef USE_EIGEN
    esu3_t ein=SU3_ECAST(in);
    SU3_ECAST(out)=SelfAdjointEigenSolver<esu3_t>(ein*ein.adjoint()).operatorInverseSqrt()*ein;
#else
    CRASH("need eigen");
#endif
  }
  
  //return a single link after the overrelaxation procedure
  void su3_find_overrelaxed(su3 out,const su3 in,const su3 staple,int nov_hits)
  {
    //compute the original contribution to the action due to the given link
    su3 prod;
    unsafe_su3_prod_su3_dag(prod,in,staple);
    
    //copy in link to out
    if(out!=in) su3_copy(out,in);
    
    //iterate over overrelax hits
    for(int ihit=0;ihit<nov_hits;ihit++)
      //scan all the three possible subgroups
      for(int isub_gr=0;isub_gr<NCOL;isub_gr++)
	{
	  //take the part of the su3 matrix
	  double r0,r1,r2,r3;
	  su2_part_of_su3(r0,r1,r2,r3,prod,isub_gr);
	  
	  //build the changing matrix
	  double x0,x1,x2,x3;
	  su2_get_overrelaxing(x0,x1,x2,x3,r0,r1,r2,r3);
	  
	  //change the link and optate the product
	  su2_prodassign_su3(x0,x1,x2,x3,isub_gr,prod);
	  su2_prodassign_su3(x0,x1,x2,x3,isub_gr,out);
	}
  }
  
  //overrelax the link
  void su3_overrelax(su3 out,const su3 in,const double w,const double *coeff,int ord)
  {
    su3 t[ord];
    
    //subtract 1 from in
    su3 f;
    su3_summ_real(f,in,-1);
    
    //ord 0
    su3_put_to_id(out);       //output init
    
    //ord 1
    su3_copy(t[1],f);
    su3_summ_the_prod_double(out,t[1],coeff[1]);
    
    //ord 2-ord
    for(int iord=2;iord<ord;iord++)
      {
	unsafe_su3_prod_su3(t[iord],t[iord-1],f);
	su3_summ_the_prod_double(out,t[iord],coeff[iord]);
      }
    
    //unitarize
    su3_unitarize_orthonormalizing(out,out);
  }
  
  //exact exponential of i times the *****passed hermitian matrix Q*****
  //algorithm taken from hep­lat/0311018
  //the stored f are relative to c0
  CUDA_HOST_AND_DEVICE void hermitian_exact_i_exponentiate_ingredients(hermitian_exp_ingredients &out,
								       const su3& Q)
  {
    //copy Q
    su3_copy(out.Q,Q);
    
  }
  
  //return a cooled copy of the passed link
  void su3_find_cooled_eo_conf(su3 u,eo_ptr<quad_su3> eo_conf,int par,int ieo,int mu)
    {
      CRASH("reimplement");//link failing
    //   //compute the staple
    // su3 staple;
    // compute_point_summed_squared_staples_eo_conf_single_dir(staple,eo_conf,loclx_of_loceo[par][ieo],mu);
    
    // //find the link that maximize the plaquette
    // su3_unitarize_maximal_trace_projecting(u,staple);
  }
  
  inline void su3_find_cooled_lx_conf(su3 u,quad_su3 *lx_conf,int ivol,int mu)
  {
    CRASH("reimplement");
    // //compute the staple
    // su3 staple;
    // compute_point_summed_squared_staples_lx_conf_single_dir(staple,lx_conf,ivol,mu);
    
    // //find the link that maximize the plaquette
    // su3_unitarize_maximal_trace_projecting(u,staple);
  }
}
