#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <string.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"

#include "complex.hpp"
#include "new_types_definitions.hpp"
#include "float_128.hpp"
#include "su3.hpp"

namespace nissa
{
  //exact exponential of i times the passed anti-hermitian matrix Q
  //algorithm taken from hep­lat/0311018
  //the stored f are relative to c0
  void anti_hermitian_exact_i_exponentiate_ingredients(anti_hermitian_exp_ingredients &out,su3 Q)
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
  
  //build the exponential from the ingredients
  void safe_anti_hermitian_exact_i_exponentiate(su3 out,anti_hermitian_exp_ingredients &ing)
  {
    //compute out according to (eq. 13)
    su3_put_to_diag(out,ing.f[0]);
    su3_summ_the_prod_complex(out,ing.Q,ing.f[1]);
    su3_summ_the_prod_complex(out,ing.Q2,ing.f[2]);
  }
  void safe_anti_hermitian_exact_i_exponentiate(su3 out,su3 Q)
  {
    anti_hermitian_exp_ingredients ing;
    anti_hermitian_exact_i_exponentiate_ingredients(ing,Q);
    safe_anti_hermitian_exact_i_exponentiate(out,ing);
  }
  
  //can be used for an hermitian matrix
  void safe_hermitian_exact_exponentiate(su3 out,su3 in)
  {
    su3 Q;
    su3_prod_idouble(Q,in,-1);
    
    safe_anti_hermitian_exact_i_exponentiate(out,Q);
  }
  
  //return sqrt(|U*U^+-1|)
  double su3_get_non_unitariness(su3 u)
  {
    su3 zero;
    su3_put_to_id(zero);
    su3_subt_the_prod_su3_dag(zero,u,u);
    
    return sqrt(real_part_of_trace_su3_prod_su3_dag(zero,zero));
  }
}
