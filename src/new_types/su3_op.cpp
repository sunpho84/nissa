#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "su3_op.hpp"
#include "geometry/geometry_eo.hpp"

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
	      crash("%lg",rotating_norm);
	  }
      }
    while(not converged);
    
    su3_copy(out,U);
  }
  
  //unitarize returning (VV^\dagger)^(-1/2)*V hep-lat/0610092
  void su3_unitarize_with_sqrt(su3 out,const su3 in)
  {
// #ifdef USE_EIGEN
//     esu3_t ein=SU3_ECAST(in);
//     SU3_ECAST(out)=SelfAdjointEigenSolver<esu3_t>(ein*ein.adjoint()).operatorInverseSqrt()*ein;
// #else
    crash("need eigen");
// #endif
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
    
    //compute the real part of the determinant (eq. 14)
    double c0=out.c0=su3_real_det(Q);
    
    //takes the square of Q
    unsafe_su3_prod_su3(out.Q2,Q,Q);
    
    //takes 1/2 of the real part of the trace of Q2 (eq. 15)
    double c1=out.c1=su3_real_trace(out.Q2)/2;
    //compute c0_max (eq. 17)
    double c0_max=2*pow(c1/3,1.5);
    
    //consider the case in which c1<4*10^-3 apart, as done in MILC
    if(c1<4e-3)
      {
	out.f[0][RE]=1-c0*c0/720;
	out.f[0][IM]=-c0*(1-c1*(1-c1/42)/20)/6;
	out.f[1][RE]=c0*(1-c1*(1-3*c1/112)/15)/24;
	out.f[1][IM]=1-c1*(1-c1*(1-c1/42)/20)/6-c0*c0/5040;
	out.f[2][RE]=0.5*(-1+c1*(1-c1*(1-c1/56)/30)/12+c0*c0/20160);
	out.f[2][IM]=0.5*(c0*(1-c1*(1-c1/48)/21)/60);
      }
    else
      {
	//take c0 module and write separately its sign (see note after eq. 34)
	out.sign=0;
	if(c0<0)
	  {
	    out.sign=1;
	    c0=-c0;
	  }
	
	//check rounding error
	double eps=(c0_max-c0)/c0_max;
	
	//(eqs. 23-24)
	double theta;
	if(eps<0) theta=out.theta=0; //only possible as an effect of rounding error when c0/c0_max=1
	else
	  if(eps<1e-3) theta=out.theta=sqrt(2*eps)*(1+(1.0/12+(3.0/160+(5.0/896+(35.0/18432+63.0/90112*eps)*eps)*eps)*eps)*eps);
	  else theta=out.theta=acos(c0/c0_max);
	double u=out.u=sqrt(c1/3)*cos(theta/3);
	double w=out.w=sqrt(c1)*sin(theta/3);
	
	//auxiliary variables for the computation of h0, h1, h2
	double u2=u*u,w2=w*w,u2mw2=u2-w2,w2p3u2=w2+3*u2,w2m3u2=w2-3*u2;
	double cu=out.cu=cos(u),c2u=out.c2u=cos(2*u);
	double su=out.su=sin(u),s2u=out.s2u=sin(2*u);
	double cw=out.cw=cos(w);
	
	//compute xi function defined after (eq. 33)
	double xi0w;
	if(fabs(w)<0.05)
	  {
	    double temp0=w*w,temp1=1-temp0/42,temp2=1.0-temp0/20*temp1;
	    xi0w=1-temp0/6*temp2;
	  }
	else xi0w=sin(w)/w;
	out.xi0w=xi0w;
	
	//computation of h0, h1, h2 (eqs. 30-32)
	complex h0={
	  u2mw2*c2u+ //(u2-w2)*cos(2u)
	  cu*8*u2*cw+ //cos(u)*8*u2*cos(w)
	  2*su*u*w2p3u2*xi0w, //sin(u)*2*mu*(3*u2+w2)*xi0(w)
	  u2mw2*s2u+ //(u2-w2)*sin(2u)
	  -su*8*u2*cw+ //-sin(u)*8*u2*cos(w)
	  cu*2*u*w2p3u2*xi0w}; //cos(u)*2*u*(3*u2+w2)*xi0(w)
	complex h1={
	  2*u*c2u+ //2*u*cos(2u)
	  -cu*2*u*cw+ //cos(u)*2*u*cos(w)
	  -su*w2m3u2*xi0w, //sin(u)*(u2-3*w2)*xi0(w)
	  2*u*s2u+ //2*u*sin(2u)
	  su*2*u*cos(w)+ //sin(u)*2*u*cos(w)
	  -cu*w2m3u2*xi0w};//cos(u)*(3*u2-w2)*xi0(w)
	complex h2={
	  c2u+ //cos(2u)
	  -cu*cw+ //-cos(u)*cos(w)
	  -3*su*u*xi0w, //-3*sin(u)*u*xi0(w)
	  s2u+ //sin(2u)
	  su*cw+ //sin(w)*cos(w)
	  -cu*3*u*xi0w};//-cos(u)*3*u*xi0(w)
	
	//build f (eq. 29)
	double fact=1/(9*u*u-w*w);
	complex_prod_double(out.f[0],h0,fact);
	complex_prod_double(out.f[1],h1,fact);
	complex_prod_double(out.f[2],h2,fact);
	
	//change sign to f according to (eq. 34)
	if(out.sign!=0)
	  {
	    out.f[0][IM]*=-1;
	    out.f[1][RE]*=-1;
	    out.f[2][IM]*=-1;
	  }
      }
  }
}
