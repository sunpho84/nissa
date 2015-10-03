#ifndef _SU3_HPP
#define _SU3_HPP

#include <math.h>
#include <string.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/random.hpp"

#include "complex.hpp"
#include "float_128.hpp"
#include "new_types_definitions.hpp"
#include "spin.hpp"
#include "su3.hpp"

#include "operations/su3_paths/squared_staples.hpp"
#include "routines/ios.hpp"
#include "routines/math_routines.hpp"

namespace nissa
{
  inline void color_put_to_zero(color m) {for(size_t ic=0;ic<NCOL;ic++) complex_put_to_zero(m[ic]);}
  inline void su3_put_to_zero(su3 m) {for(size_t ic=0;ic<NCOL;ic++) color_put_to_zero(m[ic]);}
  inline void as2t_su3_put_to_zero(as2t_su3 m) {for(size_t i=0;i<sizeof(as2t_su3)/sizeof(su3);i++) su3_put_to_zero(m[i]);}
  inline void spincolor_put_to_zero(spincolor m) {for(size_t id=0;id<4;id++) color_put_to_zero(m[id]);}
  inline void colorspinspin_put_to_zero(colorspinspin m) {for(size_t ic=0;ic<NCOL;ic++) spinspin_put_to_zero(m[ic]);}
  inline void su3spinspin_put_to_zero(su3spinspin m) {for(size_t ic=0;ic<NCOL;ic++) colorspinspin_put_to_zero(m[ic]);}
  inline void su3_put_to_id(su3 m) {su3_put_to_zero(m);for(size_t ic=0;ic<NCOL;ic++) m[ic][ic][RE]=1;}
  inline void su3_put_to_diag(su3 m,color in) {su3_put_to_zero(m);for(size_t ic=0;ic<NCOL;ic++) complex_copy(m[ic][ic],in[ic]);}
  inline void su3_put_to_diag(su3 m,complex in) {su3_put_to_zero(m);for(size_t ic=0;ic<NCOL;ic++) complex_copy(m[ic][ic],in);}
  inline void su3_put_to_diag(su3 m,double in) {su3_put_to_zero(m);for(size_t ic=0;ic<NCOL;ic++) m[ic][ic][0]=in;}
  
  //////////////////////////////////////// Copy /////////////////////////////////////
  
  inline void color_copy(color b,color a) {for(size_t ic=0;ic<NCOL;ic++) complex_copy(b[ic],a[ic]);}
  inline void su3_copy(su3 b,su3 a) {for(size_t ic=0;ic<NCOL;ic++) color_copy(b[ic],a[ic]);}
  inline void quad_su3_copy(quad_su3 b,quad_su3 a) {for(size_t i=0;i<NDIM;i++) su3_copy(b[i],a[i]);}
  inline void spincolor_copy(spincolor b,spincolor a) {for(size_t i=0;i<4;i++) color_copy(b[i],a[i]);}
  inline void colorspinspin_copy(colorspinspin b,colorspinspin a) {for(size_t i=0;i<NCOL;i++) spinspin_copy(b[i],a[i]);}
  inline void su3spinspin_copy(su3spinspin b,su3spinspin a) {for(size_t i=0;i<NCOL;i++) colorspinspin_copy(b[i],a[i]);}
  
  //////////////////// Switch directions so to agree to ILDG ordering ////////////////////
  
  inline void quad_su3_nissa_to_ildg_reord(quad_su3 out,quad_su3 in)
  {
    quad_su3 buff;
    quad_su3_copy(buff,in);
    
    for(int mu=0;mu<NDIM;mu++) su3_copy(out[(mu+NDIM-1)%NDIM],buff[mu]);
  }
  
  inline void quad_su3_ildg_to_nissa_reord(quad_su3 out,quad_su3 in)
  {
    quad_su3 buff;
    quad_su3_copy(buff,in);
    
    for(int mu=0;mu<NDIM;mu++) su3_copy(out[(mu+1)%NDIM],buff[mu]);
  }
  
  ////////////////////////////////// Operations between colors //////////////////////////
  
  //just print a color
  inline void color_print(color c)
  {
    for(size_t ic=0;ic<NCOL;ic++) printf("%+016.16lg,%+016.16lg\t",c[ic][RE],c[ic][IM]);
    master_printf("\n");
  }
  
  //summ two colors
  inline void color_summ(color a,color b,color c) {for(size_t ic=0;ic<NCOL;ic++) complex_summ(a[ic],b[ic],c[ic]);}
  inline void color_isumm(color a,color b,color c) {for(size_t ic=0;ic<NCOL;ic++) complex_isumm(a[ic],b[ic],c[ic]);}
  inline void color_isubt(color a,color b,color c) {for(size_t ic=0;ic<NCOL;ic++) complex_isubt(a[ic],b[ic],c[ic]);}
  
  inline void color_subt(color a,color b,color c) {for(size_t ic=0;ic<NCOL;ic++) complex_subt(a[ic],b[ic],c[ic]);}
  inline void color_summassign(color a,color b) {color_summ(a,a,b);}
  inline void color_subtassign(color a,color b) {color_subt(a,a,b);}
  
  inline void color_isummassign(color a,color b) {color_isumm(a,a,b);}
  inline void color_isubtassign(color a,color b) {color_isubt(a,a,b);}
  
  inline void color_prod_double(color a,color b,double c) {for(size_t ic=0;ic<NCOL;ic++) complex_prod_double(a[ic],b[ic],c);}
  inline void color_prod_idouble(color a,color b,double c) {for(size_t ic=0;ic<NCOL;ic++) complex_prod_idouble(a[ic],b[ic],c);}
  
  inline void color_summ_the_prod_double(color a,color b,double c) {for(size_t ic=0;ic<NCOL;ic++) complex_summ_the_prod_double(a[ic],b[ic],c);}
  inline void color_summ_the_prod_idouble(color a,color b,double c) {for(size_t ic=0;ic<NCOL;ic++) complex_summ_the_prod_idouble(a[ic],b[ic],c);}
  
  inline void color_summ_the_prod_complex(color a,color b,complex c) {for(size_t ic=0;ic<NCOL;ic++) complex_summ_the_prod(a[ic],b[ic],c);}
  
  inline void color_linear_comb(color a,color b,double cb,color c,double cc) {for(size_t ic=0;ic<NCOL;ic++) complex_linear_comb(a[ic],b[ic],cb,c[ic],cc);}
  
  inline void color_scalar_prod(complex out,color a,color b)
  {
    unsafe_complex_conj2_prod(out,a[0],b[0]);
    for(size_t ic=1;ic<NCOL;ic++) complex_summ_the_conj2_prod(out,a[ic],b[ic]);
  }
  inline double color_norm2(color c)
  {double out=complex_norm2(c[0]);for(size_t ic=1;ic<NCOL;ic++) out+=complex_norm2(c[ic]);return out;}
    
  /////////////////////////////// Generate an hermitean matrix ///////////////////////
  
  //Taken from M.D'Elia
#if NCOL == 3
  inline void herm_put_to_gauss(su3 H,rnd_gen *gen,double sigma)
  {
    const double one_by_sqrt3=0.577350269189626;
    const double two_by_sqrt3=1.15470053837925;
    
    double r[8];
    for(size_t ir=0;ir<4;ir++)
      {
	complex rc,ave={0,0};
	rnd_get_gauss_complex(rc,gen,ave,sigma);
	r[ir*2+0]=rc[0];
	r[ir*2+1]=rc[1];
      }
    
    //real part of diagonal elements
    H[0][0][0]= r[2]+one_by_sqrt3*r[7];
    H[1][1][0]=-r[2]+one_by_sqrt3*r[7];
    H[2][2][0]=     -two_by_sqrt3*r[7];
    
    //put immaginary part of diagonal elements to 0
    H[0][0][1]=H[1][1][1]=H[2][2][1]=0;
    
    //remaining
    H[0][1][0]=H[1][0][0]=r[0];
    H[0][1][1]=-(H[1][0][1]=r[1]);
    H[0][2][0]=H[2][0][0]=r[3];
    H[0][2][1]=-(H[2][0][1]=r[4]);
    H[1][2][0]=H[2][1][0]=r[5];
    H[1][2][1]=-(H[2][1][1]=r[6]);
  }
#else
#endif
  
  // A gauss vector has complex components z which are gaussian distributed
  // with <z~ z> = sigma
  inline void color_put_to_gauss(color H,rnd_gen *gen,double sigma)
  {
    complex ave={0,0};
    for(size_t ic=0;ic<NCOL;ic++) rnd_get_gauss_complex(H[ic],gen,ave,sigma);
  }
  
  //////////////////////////////////// Color and complex //////////////////////////////
  
  inline void safe_color_prod_complex(color out,color in,complex factor)
  {for(size_t i=0;i<NCOL;i++) safe_complex_prod(((complex*)out)[i],((complex*)in)[i],factor);}
  inline void unsafe_color_prod_complex(color out,color in,complex factor)
  {for(size_t i=0;i<NCOL;i++) unsafe_complex_prod(((complex*)out)[i],((complex*)in)[i],factor);}
  
  inline void safe_color_prod_complex_conj(color out,color in,complex factor)
  {for(size_t i=0;i<NCOL;i++) safe_complex_conj2_prod(((complex*)out)[i],((complex*)in)[i],factor);}
  inline void unsafe_color_prod_complex_conj(color out,color in,complex factor)
  {for(size_t i=0;i<NCOL;i++) unsafe_complex_conj2_prod(((complex*)out)[i],((complex*)in)[i],factor);}
  
  ////////////////////////////////// Operations between su3 //////////////////////////
  
  //comput third row according to other 2
  inline void su3_build_third_row(su3 o)
  {
    CRASH_IF_NOT_3COL();
    
    unsafe_complex_conj_conj_prod(o[2][0],o[0][1],o[1][2]);
    complex_subt_the_conj_conj_prod(o[2][0],o[0][2],o[1][1]);
    unsafe_complex_conj_conj_prod(o[2][1],o[0][2],o[1][0]);
    complex_subt_the_conj_conj_prod(o[2][1],o[0][0],o[1][2]);
    unsafe_complex_conj_conj_prod(o[2][2],o[0][0],o[1][1]);
    complex_subt_the_conj_conj_prod(o[2][2],o[0][1],o[1][0]);
  }
  
  //just print an su3 matrix
  inline void su3_print(su3 U)
  {
    for(size_t ic1=0;ic1<NCOL;ic1++)
      {
	for(size_t ic2=0;ic2<NCOL;ic2++) printf("%+016.16lg,%+016.16lg\t",U[ic1][ic2][RE],U[ic1][ic2][IM]);
	printf("\n");
      }
    printf("\n");
  }
  
  //return the trace of an su3 matrix
  inline void su3_trace(complex tr,su3 m)
  {
    complex_copy(tr,m[0][0]);
    for(size_t ic=1;ic<NCOL;ic++) complex_summassign(tr,m[ic][ic]);
  }
  
  //return only the real part of an su3 matrix
  inline double su3_real_trace(su3 m)
  {
    double a=m[0][0][RE];
    for(size_t ic=1;ic<NCOL;ic++) a+=m[ic][ic][RE];
    return a;
  }
  
  //take projection of the su2 matrix over an su3 matrix
  //return the inverse modulo of the part parallel in the original matrix
  inline double su2_part_of_su3(double &A,double &B,double &C,double &D,su3 in,int isub_gr)
  {
    //take indices of the subgroup
    int a=su3_sub_gr_indices[isub_gr][0];
    int b=su3_sub_gr_indices[isub_gr][1];
    
    //extract part parallel to sigmas
    A=in[a][a][RE]+in[b][b][RE];
    B=in[a][b][IM]+in[b][a][IM];
    C=in[a][b][RE]-in[b][a][RE];
    D=in[a][a][IM]-in[b][b][IM];
    
    //normalize
    double N=sqrt(A*A+B*B+C*C+D*D);
    if(fabs(N)==0) N=A=1;
    else
      {
	N=1/N;
	A*=N;
	B*=N;
	C*=N;
	D*=N;
      }
    
    return N;
  }
  
  //return the same in a matrix
  inline double su2_part_of_su3(su2 out,su3 in,int isub_gr)
  {
    double A,B,C,D;
    double N=su2_part_of_su3(A,B,C,D,in,isub_gr);
    
    out[0][0][RE]=A;
    out[0][0][IM]=-D;
    out[0][1][RE]=-C;
    out[0][1][IM]=-B;
    out[1][0][RE]=C;
    out[1][0][IM]=-B;
    out[1][1][RE]=A;
    out[1][1][IM]=D;
    
    return N;
  }
  
  //multiply an su2 matrix and an su3 and assign to last
  inline void su2_prodassign_su3(su2 mod,int isub_gr,su3 in)
  {
    int ic1=su3_sub_gr_indices[isub_gr][0];
    int ic2=su3_sub_gr_indices[isub_gr][1];
    
    //create the two new rows of the matrix, column by column
    for(size_t ic=0;ic<NCOL;ic++)
      {
	complex row[2];
	for(int irow=0;irow<2;irow++)
	  {
	    safe_complex_prod    (row[irow],mod[irow][0],in[ic1][ic]);
	    complex_summ_the_prod(row[irow],mod[irow][1],in[ic2][ic]);
	  }
	
	//change the two lines in the matrix
	complex_copy(in[ic1][ic],row[0]);
	complex_copy(in[ic2][ic],row[1]);
      }
  }
  
  //in the form A+B*i*sigma1+...
  inline void su2_prodassign_su3(double A,double B,double C,double D,int isub_gr,su3 in)
  {
    su2 mod={{{A,D},{C,B}},{{-C,B},{A,-D}}};
    su2_prodassign_su3(mod,isub_gr,in);
  }
  
  //summ the trace to the input
  inline void su3_summ_the_trace(complex tr,su3 m)
  {for(size_t ic=0;ic<NCOL;ic++) complex_summ(tr,tr,m[ic][ic]);}
  
  //return the anti-hermitian traceless part of an su3 matrix
  inline void unsafe_su3_traceless_anti_hermitian_part(su3 out,su3 in)
  {
    double trace_im_third=0;
    for(size_t ic=0;ic<NCOL;ic++) trace_im_third+=in[ic][ic][IM];
    trace_im_third/=NCOL;
    
    for(size_t ic=0;ic<NCOL;ic++)
      {
	//real part of diagonal: 0
	out[ic][ic][RE]=0;
	//imag part of diagonal: subtract the trace
	out[ic][ic][IM]=in[ic][ic][IM]-trace_im_third;
	
	for(size_t jc=0;jc<ic;jc++)
	  {
	    //out-of-diag real part
	    out[ic][jc][RE]=-(out[jc][ic][RE]=(in[jc][ic][RE]-in[ic][jc][RE])/2);
	    //out-of-diag imag part
	    out[ic][jc][IM]=+(out[jc][ic][IM]=(in[jc][ic][IM]+in[ic][jc][IM])/2);
	  }
      }
  }
  
  //return the hermitian traceless part of an su3 matrix
  inline void unsafe_su3_traceless_hermitian_part(su3 out,su3 in)
  {
    double trace_re_third=0;
    for(size_t ic=0;ic<NCOL;ic++) trace_re_third+=in[ic][ic][RE];
    trace_re_third/=NCOL;
    
    for(size_t ic=0;ic<NCOL;ic++)
      {
	//imag part of diagonal: 0
	out[ic][ic][IM]=0;
	//real part of diagonal: subtract the trace
	out[ic][ic][RE]=in[ic][ic][RE]-trace_re_third;
	
	for(size_t jc=0;jc<ic;jc++)
	  {
	    //out-of-diag real part
	    out[ic][jc][RE]=+(out[jc][ic][RE]=(in[jc][ic][RE]+in[ic][jc][RE])/2);
	    //out-of-diag imag part
	    out[ic][jc][IM]=-(out[jc][ic][IM]=(in[jc][ic][IM]-in[ic][jc][IM])/2);
	  }
      }
  }
  
  //calculate the determinant of an su3 matrix
  inline void su3_det(complex d,su3 U)
  {
#if NCOL == 3
    complex a;
    
    unsafe_complex_prod(a,U[1][1],U[2][2]);
    complex_subt_the_prod(a,U[1][2],U[2][1]);
    unsafe_complex_prod(d,U[0][0],a);
    
    unsafe_complex_prod(a,U[1][2],U[2][0]);
    complex_subt_the_prod(a,U[1][0],U[2][2]);
    complex_summ_the_prod(d,U[0][1],a);
    
    unsafe_complex_prod(a,U[1][0],U[2][1]);
    complex_subt_the_prod(a,U[1][1],U[2][0]);
    complex_summ_the_prod(d,U[0][2],a);
#else
    matrix_determinant(d,(complex*)U,NCOL);
#endif
  }
  
  //calculate the real part of the determinant of an su3 matrix
  inline double su3_real_det(su3 u)
  {
#if NCOL == 3
    return
      u[0][2][IM]*(u[1][1][RE]*u[2][0][IM]+u[1][1][IM]*u[2][0][RE]-u[1][0][RE]*u[2][1][IM]-u[1][0][IM]*u[2][1][RE])+
      u[0][2][RE]*(u[1][1][IM]*u[2][0][IM]-u[1][1][RE]*u[2][0][RE]-u[1][0][IM]*u[2][1][IM]+u[1][0][RE]*u[2][1][RE])+
      u[0][1][IM]*(u[1][0][RE]*u[2][2][IM]-u[1][2][RE]*u[2][0][IM]-u[1][2][IM]*u[2][0][RE]+u[1][0][IM]*u[2][2][RE])+
      u[0][1][RE]*(u[1][0][IM]*u[2][2][IM]-u[1][2][IM]*u[2][0][IM]+u[1][2][RE]*u[2][0][RE]-u[1][0][RE]*u[2][2][RE])+
      u[0][0][IM]*(u[1][2][RE]*u[2][1][IM]+u[1][2][IM]*u[2][1][RE]-u[1][1][RE]*u[2][2][IM]-u[1][1][IM]*u[2][2][RE])+
      u[0][0][RE]*(u[1][2][IM]*u[2][1][IM]-u[1][2][RE]*u[2][1][RE]-u[1][1][IM]*u[2][2][IM]+u[1][1][RE]*u[2][2][RE]);
#else
    complex d;
    su3_det(d,u);
    return d[RE];
#endif
  }
  
  //return the hemitian su3 matrix
  inline void unsafe_su3_hermitian(su3 out,su3 in)
  {
    for(size_t ic_in=0;ic_in<NCOL;ic_in++)
      for(size_t ic_out=0;ic_out<NCOL;ic_out++)
	complex_conj(out[ic_in][ic_out],in[ic_out][ic_in]);
  }
  inline void safe_su3_hermitian(su3 out,su3 in)
  {
    su3 tmp;
    unsafe_su3_hermitian(tmp,in);
    su3_copy(out,tmp);
  }
  
  //return the transposed su3 matrix
  inline void unsafe_su3_transpose(su3 out,su3 in)
  {
    for(size_t ic_in=0;ic_in<NCOL;ic_in++)
      for(size_t ic_out=0;ic_out<NCOL;ic_out++)
	complex_copy(out[ic_in][ic_out],in[ic_out][ic_in]);
  }
  inline void safe_su3_transpose(su3 out,su3 in)
  {
    su3 tmp;
    unsafe_su3_transpose(tmp,in);
    su3_copy(out,tmp);
  }
  
  //summ two su3 matrixes
  inline void su3_summ(su3 a,su3 b,su3 c) {for(size_t ic=0;ic<NCOL;ic++) color_summ(a[ic],b[ic],c[ic]);}
  inline void unsafe_su3_summ_su3_dag(su3 a,su3 b,su3 c) {for(size_t i=0;i<NCOL;i++) for(size_t j=0;j<NCOL;j++) complex_summ_conj2(a[i][j],b[i][j],c[j][i]);}
  inline void su3_summassign_su3_dag(su3 a,su3 b) {unsafe_su3_summ_su3_dag(a,a,b);}
  inline void su3_summassign(su3 a,su3 b){su3_summ(a,a,b);}
  inline void su3_summ_real(su3 a,su3 b,double c) {su3_copy(a,b);for(size_t i=0;i<NCOL;i++) a[i][i][0]=b[i][i][0]+c;}
  inline void su3_subt(su3 a,su3 b,su3 c) {for(size_t ic=0;ic<NCOL;ic++) color_subt(a[ic],b[ic],c[ic]);}
  inline void su3_subtassign(su3 a,su3 b) {su3_subt(a,a,b);}
  inline void su3_subt_complex(su3 a,su3 b,complex c) {su3_copy(a,b);for(size_t i=0;i<NCOL;i++) complex_subt(a[i][i],b[i][i],c);}
  inline void unsafe_su3_subt_su3_dag(su3 a,su3 b,su3 c) {for(size_t i=0;i<NCOL;i++) for(size_t j=0;j<NCOL;j++) complex_subt_conj2(a[i][j],b[i][j],c[j][i]);}
  inline void su3_subtassign_su3_dag(su3 a,su3 b) {unsafe_su3_subt_su3_dag(a,a,b);}
  
  //Product of two su3 matrixes
  inline void unsafe_su3_prod_su3(su3 a,su3 b,su3 c,const size_t nr_max=NCOL)
  {
    for(size_t ir_out=0;ir_out<nr_max;ir_out++)
      for(size_t ic_out=0;ic_out<NCOL;ic_out++)
	{
	  unsafe_complex_prod(a[ir_out][ic_out],b[ir_out][0],c[0][ic_out]);
	  for(size_t itemp=1;itemp<NCOL;itemp++)
	    complex_summ_the_prod(a[ir_out][ic_out],b[ir_out][itemp],c[itemp][ic_out]);
	}
  }
  inline void safe_su3_prod_su3(su3 a,su3 b,su3 c) {su3 d;unsafe_su3_prod_su3(d,b,c);su3_copy(a,d);}
  inline void su3_prodassign_su3(su3 a,su3 b) {safe_su3_prod_su3(a,a,b);}
  inline void su3_summ_the_prod_su3(su3 a,su3 b,su3 c)
  {
    for(size_t ir_out=0;ir_out<NCOL;ir_out++)
      for(size_t ic_out=0;ic_out<NCOL;ic_out++)
	for(size_t itemp=0;itemp<NCOL;itemp++)
	  complex_summ_the_prod(a[ir_out][ic_out],b[ir_out][itemp],c[itemp][ic_out]);
  }
  inline void su3_summ_the_prod_su3_dag(su3 a,su3 b,su3 c)
  {
    for(size_t ir_out=0;ir_out<NCOL;ir_out++)
      for(size_t ic_out=0;ic_out<NCOL;ic_out++)
	for(size_t itemp=0;itemp<NCOL;itemp++)
	  complex_summ_the_conj2_prod(a[ir_out][ic_out],b[ir_out][itemp],c[ic_out][itemp]);
  }
  inline void su3_summ_the_dag_prod_su3(su3 a,su3 b,su3 c)
  {
    for(size_t ir_out=0;ir_out<NCOL;ir_out++)
      for(size_t ic_out=0;ic_out<NCOL;ic_out++)
	for(size_t itemp=0;itemp<NCOL;itemp++)
	  complex_summ_the_conj1_prod(a[ir_out][ic_out],b[itemp][ir_out],c[itemp][ic_out]);
  }
  
  //Product of two su3 matrixes
  inline void unsafe_su3_dag_prod_su3(su3 a,su3 b,su3 c,const size_t nr_max=NCOL)
  {
    for(size_t ir_out=0;ir_out<nr_max;ir_out++)
      for(size_t ic_out=0;ic_out<NCOL;ic_out++)
	{
	  unsafe_complex_conj1_prod(a[ir_out][ic_out],b[0][ir_out],c[0][ic_out]);
	  for(size_t itemp=1;itemp<NCOL;itemp++)
	    complex_summ_the_conj1_prod(a[ir_out][ic_out],b[itemp][ir_out],c[itemp][ic_out]);
	}
  }
  inline void safe_su3_dag_prod_su3(su3 a,su3 b,su3 c) {su3 d;unsafe_su3_dag_prod_su3(d,b,c);su3_copy(a,d);}
  inline void su3_dag_summ_the_prod_su3(su3 a,su3 b,su3 c)
  {
    for(size_t ir_out=0;ir_out<NCOL;ir_out++)
      for(size_t ic_out=0;ic_out<NCOL;ic_out++)
	for(size_t itemp=0;itemp<NCOL;itemp++)
	  complex_summ_the_conj1_prod(a[ir_out][ic_out],b[itemp][ir_out],c[itemp][ic_out]);
  }
  
  //Product of two su3 matrixes
  inline void unsafe_su3_prod_su3_dag(su3 a,su3 b,su3 c,const size_t nr_max=NCOL)
  {
    for(size_t ir_out=0;ir_out<nr_max;ir_out++)
      for(size_t ic_out=0;ic_out<NCOL;ic_out++)
	{
	  unsafe_complex_conj2_prod(a[ir_out][ic_out],b[ir_out][0],c[ic_out][0]);
	  for(size_t jc=1;jc<NCOL;jc++) complex_summ_the_conj2_prod(a[ir_out][ic_out],b[ir_out][jc],c[ic_out][jc]);
	}
  }
  inline void safe_su3_prod_su3_dag(su3 a,su3 b,su3 c) {su3 d;unsafe_su3_prod_su3_dag(d,b,c);su3_copy(a,d);}
  
  //subtract the product
  inline void su3_subt_the_prod_su3_dag(su3 a,su3 b,su3 c)
  {
    for(size_t ir_out=0;ir_out<NCOL;ir_out++)
      for(size_t ic_out=0;ic_out<NCOL;ic_out++)
	for(size_t jc=0;jc<NCOL;jc++)
	  complex_subt_the_conj2_prod(a[ir_out][ic_out],b[ir_out][jc],c[ic_out][jc]);
  }
  
  //Trace of the product of two su3 matrices
  inline double real_part_of_trace_su3_prod_su3_dag(su3 a,su3 b)
  {
    double t=0;
    
    for(size_t ic1=0;ic1<NCOL;ic1++)
      for(size_t ic2=0;ic2<NCOL;ic2++)
	t+=real_part_of_complex_scalar_prod(a[ic1][ic2],b[ic1][ic2]);
    
    return t;
  }
  
  //Trace of the product of two su3 matrices
  inline void trace_su3_prod_su3(complex t,su3 a,su3 b)
  {
    complex_put_to_zero(t);
    for(size_t ic1=0;ic1<NCOL;ic1++)
      for(size_t ic2=0;ic2<NCOL;ic2++)
	complex_summ_the_prod(t,a[ic1][ic2],b[ic2][ic1]);
  }
  
  //Product of two su3 matrixes
  inline void unsafe_su3_dag_prod_su3_dag(su3 a,su3 b,su3 c,const size_t nr_max=NCOL)
  {
    for(size_t ir_out=0;ir_out<nr_max;ir_out++)
      for(size_t ic_out=0;ic_out<NCOL;ic_out++)
	{
	  unsafe_complex_conj_conj_prod(a[ir_out][ic_out],b[0][ir_out],c[ic_out][0]);
	  for(size_t itemp=1;itemp<NCOL;itemp++)
	    complex_summ_the_conj_conj_prod(a[ir_out][ic_out],b[itemp][ir_out],c[ic_out][itemp]);
	}
  }
  inline void safe_su3_dag_prod_su3_dag(su3 a,su3 b,su3 c)
  {su3 d;unsafe_su3_dag_prod_su3_dag(d,b,c);su3_copy(a,d);}
  
  //product of an su3 matrix by a complex
  inline void unsafe_su3_prod_complex(su3 a,su3 b,complex c)
  {for(size_t ic=0;ic<NCOL;ic++) unsafe_color_prod_complex(a[ic],b[ic],c);}
  inline void safe_su3_prod_complex(su3 a,su3 b,complex c)
  {su3 d;unsafe_su3_prod_complex(d,b,c);su3_copy(a,d);}
  
  //product of an su3 matrix by a complex
  inline void unsafe_su3_prod_complex_conj(su3 a,su3 b,complex c)
  {for(size_t ic=0;ic<NCOL;ic++) unsafe_color_prod_complex_conj(a[ic],b[ic],c);}
  inline void safe_su3_prod_complex_conj(su3 a,su3 b,complex c)
  {su3 d;unsafe_su3_prod_complex_conj(d,b,c);su3_copy(a,d);}
  
  //product of an su3 matrix by a real
  inline void su3_prod_double(su3 a,su3 b,double r)
  {for(size_t ic=0;ic<NCOL;ic++) color_prod_double(a[ic],b[ic],r);}
  inline void su3_prodassign_double(su3 a,double r)
  {su3_prod_double(a,a,r);}
  
  //hermitian of su3 matrix times a real
  inline void unsafe_su3_hermitian_prod_double(su3 a,su3 b,double r)
  {
    for(size_t i=0;i<NCOL;i++)
      for(size_t j=0;j<NCOL;j++)
	{
	  a[i][j][RE]=+r*b[j][i][RE];
	  a[i][j][IM]=-r*b[j][i][IM];
	}
  }
  inline void safe_su3_hermitian_prod_double(su3 a,su3 b,double r)
  {
    for(size_t i=0;i<NCOL;i++)
      {
	a[i][i][RE]=+r*b[i][i][RE];
	a[i][i][IM]=-r*b[i][i][IM];
	for(size_t j=i+1;j<NCOL;j++)
	  {
	    double a_i_j_RE=+r*b[j][i][RE];
	    double a_i_j_IM=-r*b[j][i][IM];
	    a[j][i][RE]=+r*b[i][j][RE];
	    a[j][i][IM]=-r*b[i][j][IM];
	    a[i][j][RE]=a_i_j_RE;
	    a[i][j][IM]=a_i_j_IM;
	  }
      }
  }
  inline void su3_summ_the_hermitian_prod_double(su3 a,su3 b,double r)
  {
    for(size_t i=0;i<NCOL;i++)
      for(size_t j=0;j<NCOL;j++)
	{
	  a[i][j][RE]+= r*b[j][i][RE];
	  a[i][j][IM]+=-r*b[j][i][IM];
	}
  }
  
  //summ the prod of su3 with imag
  inline void su3_prod_idouble(su3 a,su3 b,double r) {for(size_t i=0;i<NCOL;i++) color_prod_idouble(a[i],b[i],r);}
  inline void su3_prodassign_idouble(su3 a,double r) {su3_prod_idouble(a,a,r);}
  inline void su3_summ_the_prod_idouble(su3 a,su3 b,double r) {for(size_t i=0;i<NCOL;i++) color_summ_the_prod_idouble(a[i],b[i],r);}
  //summ the prod of su3 with real
  inline void su3_summ_the_prod_double(su3 a,su3 b,double r) {for(size_t i=0;i<NCOL;i++) color_summ_the_prod_double(a[i],b[i],r);}
  
  //summ the prod of the dag su3 with real
  inline void su3_dag_summ_the_prod_double(su3 a,su3 b,double r)
  {
    for(size_t i=0;i<NCOL;i++)
      for(size_t j=0;j<NCOL;j++)
	{
	  a[i][j][RE]+=b[j][i][RE]*r;
	  a[i][j][IM]-=b[j][i][IM]*r;
	}
  }
  
  //combine linearly two su3 elements
  inline void su3_linear_comb(su3 a,su3 b,double cb,su3 c,double cc) {for(size_t ic=0;ic<NCOL;ic++) color_linear_comb(a[ic],b[ic],cb,c[ic],cc);}
  //summ the prod of su3 with complex
  inline void su3_summ_the_prod_complex(su3 a,su3 b,complex c) {for(size_t ic=0;ic<NCOL;ic++) color_summ_the_prod_complex(a[ic],b[ic],c);}
  
  //calculate explicitely the inverse
  inline void unsafe_su3_explicit_inverse(su3 invU,su3 U)
  {
    CRASH_IF_NOT_3COL();
    
    complex det,rec_det;
    su3_det(det,U);
    complex_reciprocal(rec_det,det);
    
    unsafe_complex_prod(invU[0][0],U[1][1],U[2][2]);
    unsafe_complex_prod(invU[1][0],U[1][2],U[2][0]);
    unsafe_complex_prod(invU[2][0],U[1][0],U[2][1]);
    
    unsafe_complex_prod(invU[0][1],U[0][2],U[2][1]);
    unsafe_complex_prod(invU[1][1],U[0][0],U[2][2]);
    unsafe_complex_prod(invU[2][1],U[0][1],U[2][0]);
    
    unsafe_complex_prod(invU[0][2],U[0][1],U[1][2]);
    unsafe_complex_prod(invU[1][2],U[0][2],U[1][0]);
    unsafe_complex_prod(invU[2][2],U[0][0],U[1][1]);
    
    
    complex_subt_the_prod(invU[0][0],U[1][2],U[2][1]);
    complex_subt_the_prod(invU[1][0],U[1][0],U[2][2]);
    complex_subt_the_prod(invU[2][0],U[1][1],U[2][0]);
    
    complex_subt_the_prod(invU[0][1],U[0][1],U[2][2]);
    complex_subt_the_prod(invU[1][1],U[0][2],U[2][0]);
    complex_subt_the_prod(invU[2][1],U[0][0],U[2][1]);
    
    complex_subt_the_prod(invU[0][2],U[0][2],U[1][1]);
    complex_subt_the_prod(invU[1][2],U[0][0],U[1][2]);
    complex_subt_the_prod(invU[2][2],U[0][1],U[1][0]);
    
    for(size_t icol1=0;icol1<NCOL;icol1++)
      for(size_t icol2=0;icol2<NCOL;icol2++)
	safe_complex_prod(invU[icol1][icol2],invU[icol1][icol2],rec_det);
  }
  inline void safe_su3_explicit_inverse(su3 invU,su3 U)
  {su3 tempU;unsafe_su3_explicit_inverse(tempU,U);su3_copy(invU,tempU);}
  
  //summ of the squared norm of the entries
  inline double su3_norm2(su3 U)
  {
    double norm2=0;
    for(size_t ic=0;ic<NCOL;ic++) norm2+=color_norm2(U[ic]);
    
    return norm2;
  }
  
  //compute the square root of y numerically
  inline void unsafe_su3_sqrt(su3 x,su3 y)
  {
    su3_put_to_id(x);
    
    double err;
    do
      {
	//update
	su3 t;
	unsafe_su3_explicit_inverse(t,x);
	safe_su3_prod_su3(t,y,t);
	su3_subt(t,t,x);
	su3_summ_the_prod_double(x,t,0.5);
	
	//compute the error
	unsafe_su3_prod_su3(t,x,x);
	su3_subt(t,y,t);
	err=sqrt(su3_norm2(t));
      }
    while(err>1.e-15);
  }
  inline void safe_su3_sqrt(su3 x,su3 y) {su3 t;unsafe_su3_sqrt(t,y);su3_copy(x,t);}
  
  //exponentiate an su3 matrix through taylor expansion
  inline void unsafe_su3_taylor_exponentiate(su3 out,su3 in,size_t order)
  {
    //1st terms
    double coef=1;
    su3 temp;
    su3_copy(temp,in);
    
    //order 0+1
    su3_put_to_id(out);
    su3_summassign(out,in);
    
    //higher orders
    for(size_t iorder=2;iorder<=order;iorder++)
      {
	safe_su3_prod_su3(temp,temp,in);
	coef/=iorder;
	su3_summ_the_prod_double(out,temp,coef);
      }
  }
  
  //unitarize by orthonormalizing the rows
  inline void su3_unitarize_orthonormalizing(su3 o,su3 i)
  {
    CRASH_IF_NOT_3COL();
    
    //compute the squared norm of row 0
    double row0_norm2=color_norm2(i[0]);
    
    //compute the scalar product between row 1 and 0
    complex row10_sc_prod;
    unsafe_complex_conj2_prod(  row10_sc_prod,i[1][0],i[0][0]);
    complex_summ_the_conj2_prod(row10_sc_prod,i[1][1],i[0][1]);
    complex_summ_the_conj2_prod(row10_sc_prod,i[1][2],i[0][2]);
    
    //orthonormalize row 1
    complex f;
    complex_prod_double(f,row10_sc_prod,1/row0_norm2);
    
    for(size_t c=0;c<NCOL;c++)
      {
	for(size_t ri=0;ri<2;ri++) o[1][c][ri]=i[1][c][ri];
	complex_subt_the_prod(o[1][c],f,i[0][c]);
      }
    
    double row0_norm=1/sqrt(row0_norm2);
    double row1_norm=1/sqrt(complex_norm2(o[1][0])+complex_norm2(o[1][1])+complex_norm2(o[1][2]));
    
    //normalize the rows
    for(size_t c=0;c<NCOL;c++)
      for(size_t ri=0;ri<2;ri++)
	{
	  o[0][c][ri]=row0_norm*i[0][c][ri];
	  o[1][c][ri]=row1_norm*o[1][c][ri];
	}
    
    //orthonormalize the third row
    unsafe_complex_conj_conj_prod(o[2][0],o[0][1],o[1][2]);
    complex_subt_the_conj_conj_prod(o[2][0],o[0][2],o[1][1]);
    unsafe_complex_conj_conj_prod(o[2][1],o[0][2],o[1][0]);
    complex_subt_the_conj_conj_prod(o[2][1],o[0][0],o[1][2]);
    unsafe_complex_conj_conj_prod(o[2][2],o[0][0],o[1][1]);
    complex_subt_the_conj_conj_prod(o[2][2],o[0][1],o[1][0]);
  }
  
  void su3_find_overrelaxed(su3 out,su3 in,su3 staple,int nov_hits);
  void su3_overrelax(su3 out,su3 in,double w,double *coeff,const int ord);
  
  //overrelax the link using approximated exponentiation
  inline void su3_overrelax(su3 out,su3 in,double w)
  {
    double coeff[5]={1,w,w*(w-1)/2,w*(w-1)*(w-2)/6,w*(w-1)*(w-2)*(w-3)/24};
    su3_overrelax(out,in,w,coeff,5);
  }
  
  //unitarize an su3 matrix by taking explicitely the inverse and averaging with it
  inline void su3_unitarize_explicitly_inverting(su3 new_link,su3 prop_link)
  {
    su3 inv,temp_link;
    double gamma,residue;
    
    su3_copy(temp_link,prop_link);
    
    do
      {
	unsafe_su3_explicit_inverse(inv,temp_link);
	gamma=sqrt(su3_norm2(inv)/su3_norm2(temp_link));
	
	//average U and U^-1^+
	residue=0;
	for(size_t icol1=0;icol1<NCOL;icol1++)
	  for(size_t icol2=0;icol2<NCOL;icol2++)
	    {
	      new_link[icol1][icol2][RE]=0.5*(temp_link[icol1][icol2][RE]*gamma+inv[icol2][icol1][RE]/gamma);
	      new_link[icol1][icol2][IM]=0.5*(temp_link[icol1][icol2][IM]*gamma-inv[icol2][icol1][IM]/gamma);
	      for(size_t ri=0;ri<2;ri++)
		{
		  double diff=new_link[icol1][icol2][ri]-temp_link[icol1][icol2][ri];
		  residue+=diff*diff;
		}
	    }
	
	su3_copy(temp_link,new_link);
	
	residue=sqrt(residue);
      }
    while(residue>1.e-15);
    
    //divide by third root of det
    complex det,fact;
    su3_det(det,new_link);
    complex_pow(fact,det,-1.0/NCOL);
    safe_su3_prod_complex(new_link,new_link,fact);
  }
  
  //perform an iteration of maximal projection trace
  //by maximing the trace over three different subgroups of su3
  inline void su3_unitarize_maximal_trace_projecting_iteration_slow(su3 U,su3 M)
  {
    //loop over the three subgroups
    for(size_t isub_gr=0;isub_gr<NCOL;isub_gr++)
      {
	//compute the product
	su3 prod;
	unsafe_su3_prod_su3_dag(prod,U,M);
	
	//take the subgroup isub_gr
	su2 sub;
	su2_part_of_su3(sub,prod,isub_gr);
	
	//modify the subgroup
	su2_prodassign_su3(sub,isub_gr,U);
      }
  }
  
  //perform an iteration of maximal projection trace
  //by maximing the trace over three different subgroups of su3
  inline double su3_unitarize_maximal_trace_projecting_iteration(su3 U,su3 M)
  {
    //compute the product
    su3 prod;
    unsafe_su3_prod_su3_dag(prod,U,M);
    
    //master_printf("isub -1, %16.16lg\n",su3_real_trace(prod));
    //loop over the three subgroups
    for(size_t isub_gr=0;isub_gr<NCOL;isub_gr++)
      {
	//take the subgroup isub_gr
	su2 sub;
	su2_part_of_su3(sub,prod,isub_gr);
	
	//modify the subgroup
	su2_prodassign_su3(sub,isub_gr,U);
	su2_prodassign_su3(sub,isub_gr,prod);
	
	//master_printf("isub %d, %16.16lg\n",isub_gr,su3_real_trace(prod));
      }
    
    return su3_real_trace(prod);
  }
  
  //perform maximal projection trace up to reaching the machine precision
  inline void su3_unitarize_maximal_trace_projecting(su3 out,su3 M)
  {
    //initialize the guess with the identity - proved to be faster than any good guess,
    //because iterations are so good
    su3 U;
    su3_put_to_id(U);
    
    //compute the "product", that means taking dag of M as U=1
    su3 prod;
    unsafe_su3_hermitian(prod,M);
    
    //traces to go out from loop
    double new_trace=su3_real_trace(prod);
    double old_trace;
    
    int iter=0;
    do
      {
	//store old trace
	old_trace=new_trace;
	
	//fix subgroup
	int isub_gr=iter%NCOL;
	
	//take the subgroup isub_gr
	su2 sub;
	su2_part_of_su3(sub,prod,isub_gr);
	
	//modify the subgroup
	su2_prodassign_su3(sub,isub_gr,U);
	
	//modify the prod and compute trace
	su2_prodassign_su3(sub,isub_gr,prod);
	new_trace=su3_real_trace(prod);
	
	iter++;
      }
    while(new_trace>old_trace);
    
    su3_copy(out,U);
  }
  
  void su3_find_heatbath(su3 out,su3 in,su3 staple,double beta,int nhb_hits,rnd_gen *gen);
  void su3_find_overrelaxed(su3 out,su3 in,su3 staple,int nov_hits);
  
  //return a cooled copy of the passed link
  inline void su3_find_cooled_eo_conf(su3 u,quad_su3 **eo_conf,int par,int ieo,int mu)
  {
    //compute the staple
    su3 staple;
    compute_point_summed_squared_staples_eo_conf_single_dir(staple,eo_conf,loclx_of_loceo[par][ieo],mu);
    
    //find the link that maximize the plaquette
    su3_unitarize_maximal_trace_projecting(u,staple);
  }
  inline void su3_find_cooled_lx_conf(su3 u,quad_su3 *lx_conf,int ivol,int mu)
  {
    //compute the staple
    su3 staple;
    compute_point_summed_squared_staples_lx_conf_single_dir(staple,lx_conf,ivol,mu);
    
    //find the link that maximize the plaquette
    su3_unitarize_maximal_trace_projecting(u,staple);
  }
  
  ////////////////////// products between su3 and color //////////////////
  
  //product of an su3 matrix by a color vector
  inline void unsafe_su3_prod_color(color a,su3 b,color c)
  {
    for(size_t c1=0;c1<NCOL;c1++)
      {
	unsafe_complex_prod(a[c1],b[c1][0],c[0]);
	for(size_t c2=1;c2<NCOL;c2++) complex_summ_the_prod(a[c1],b[c1][c2],c[c2]);
      }
  }
  
  //product of an su3 matrix by a color vector
  inline void unsafe_single_su3_prod_single_color(single_color a,single_su3 b,single_color c)
  {
    for(size_t c1=0;c1<NCOL;c1++)
      {
	unsafe_single_complex_prod(a[c1],b[c1][0],c[0]);
	for(size_t c2=1;c2<NCOL;c2++) single_complex_summ_the_prod(a[c1],b[c1][c2],c[c2]);
      }
  }
  
  //safe prod
  inline void safe_su3_prod_color(color a,su3 b,color c) {color t;unsafe_su3_prod_color(t,b,c);color_copy(a,t);}
  
  //summ
  inline void su3_summ_the_prod_color(color a,su3 b,color c)
  {for(size_t c1=0;c1<NCOL;c1++) for(size_t c2=0;c2<NCOL;c2++) complex_summ_the_prod(a[c1],b[c1][c2],c[c2]);}
  inline void single_su3_summ_the_prod_single_color(single_color a,single_su3 b,single_color c)
  {for(size_t c1=0;c1<NCOL;c1++) for(size_t c2=0;c2<NCOL;c2++) single_complex_summ_the_prod(a[c1],b[c1][c2],c[c2]);}
  //subt
  inline void su3_subt_the_prod_color(color a,su3 b,color c)
  {for(size_t c1=0;c1<NCOL;c1++) for(size_t c2=0;c2<NCOL;c2++) complex_subt_the_prod(a[c1],b[c1][c2],c[c2]);}
  
  //dag prod
  inline void unsafe_su3_dag_prod_color(color a,su3 b,color c)
  {
    for(size_t c1=0;c1<NCOL;c1++)
      {
	unsafe_complex_conj1_prod(a[c1],b[0][c1],c[0]);
	for(size_t c2=1;c2<NCOL;c2++) complex_summ_the_conj1_prod(a[c1],b[c2][c1],c[c2]);
      }
  }
  
  //safe dag prod
  inline void safe_su3_dag_prod_color(color a,su3 b,color c)
  {color t;unsafe_su3_dag_prod_color(t,b,c);color_copy(a,t);}
  
  //summ dag
  inline void su3_dag_summ_the_prod_color(color a,su3 b,color c)
  {for(size_t c1=0;c1<NCOL;c1++) for(size_t c2=0;c2<NCOL;c2++) complex_summ_the_conj1_prod(a[c1],b[c2][c1],c[c2]);}
  inline void single_su3_dag_summ_the_prod_single_color(single_color a,single_su3 b,single_color c)
  {for(size_t c1=0;c1<NCOL;c1++) for(size_t c2=0;c2<NCOL;c2++)single_complex_summ_the_conj1_prod(a[c1],b[c2][c1],c[c2]);}
  
  //subt dag
  inline void su3_dag_subt_the_prod_color(color a,su3 b,color c)
  {for(size_t c1=0;c1<NCOL;c1++) for(size_t c2=0;c2<NCOL;c2++) complex_subt_the_conj1_prod(a[c1],b[c2][c1],c[c2]);}
  inline void single_su3_dag_subt_the_prod_single_color(single_color a,single_su3 b,single_color c)
  {for(size_t c1=0;c1<NCOL;c1++) for(size_t c2=0;c2<NCOL;c2++) single_complex_subt_the_conj1_prod(a[c1],b[c2][c1],c[c2]);}
  
  //////////////////////////////////////////// color prod su3 ///////////////////////////////////////////
  
  //prod
  inline void unsafe_color_prod_su3(color a,color b,su3 c)
  {
    for(size_t c1=0;c1<NCOL;c1++)
      {
	unsafe_complex_prod(a[c1],b[0],c[0][c1]);
	for(size_t c2=1;c2<NCOL;c2++) complex_summ_the_prod(a[c1],b[c2],c[c2][c1]);
      }
  }
  inline void safe_color_prod_su3(color a,color b,su3 c)
  {color t;unsafe_color_prod_su3(t,b,c);color_copy(a,t);}
  
  //dag
  inline void unsafe_color_prod_su3_dag(color a,color b,su3 c)
  {
    for(size_t c1=0;c1<NCOL;c1++)
      {
	unsafe_complex_conj2_prod(a[c1],b[0],c[c1][0]);
	for(size_t c2=1;c2<NCOL;c2++) complex_summ_the_conj2_prod(a[c1],b[c2],c[c1][c2]);
      }
  }
  
  ////////////////////////////////// Operations between spincolor ////////////////////////
  
  //just print a spincolor
  inline void spincolor_print(spincolor c)
  {
    for(size_t id=0;id<4;id++)
      {
	for(size_t ic=0;ic<NCOL;ic++) printf("%+016.16le,%+016.16le\t",c[id][ic][0],c[id][ic][1]);
	master_printf("\n");
      }
    master_printf("\n");
  }
  
  //summ two spincolors
  inline void spincolor_summ(spincolor a,spincolor b,spincolor c) {for(size_t i=0;i<4;i++) color_summ(a[i],b[i],c[i]);}
  inline void spincolor_summassign(spincolor a,spincolor b) {spincolor_summ(a,a,b);}
  
  //subtract two spincolors
  inline void spincolor_subt(spincolor a,spincolor b,spincolor c) {for(size_t i=0;i<4;i++) color_subt(a[i],b[i],c[i]);}
  inline void spincolor_subtassign(spincolor a,spincolor b) {spincolor_subt(a,a,b);}
  
  //spincolor*real
  inline void spincolor_prod_double(spincolor a,spincolor b,double factor)
  {for(size_t i=0;i<4;i++) color_prod_double(a[i],b[i],factor);}
  
  //spincolor*complex
  inline void unsafe_spincolor_prod_complex(spincolor a,spincolor b,complex factor)
  {for(size_t i=0;i<4;i++) unsafe_color_prod_complex(a[i],b[i],factor);}
  inline void safe_spincolor_prod_complex(spincolor a,spincolor b,complex factor)
  {spincolor c;spincolor_copy(c,b);unsafe_spincolor_prod_complex(a,c,factor);}
  
  //spincolor+spincolor*complex
  inline void spincolor_summ_the_prod_complex(spincolor a,spincolor b,complex factor)
  {for(size_t i=0;i<4;i++) color_summ_the_prod_complex(a[i],b[i],factor);}
  
  //spincolor+spincolor*idouble
  inline void spincolor_summ_the_prod_idouble(spincolor a,spincolor b,double c) {for(size_t i=0;i<4;i++) color_summ_the_prod_idouble(a[i],b[i],c);}
  
  //dirac*spincolor
  inline void unsafe_dirac_prod_spincolor(spincolor out,dirac_matr *m,spincolor in)
  {for(size_t id1=0;id1<4;id1++) unsafe_color_prod_complex(out[id1],in[m->pos[id1]],m->entr[id1]);}
  
  //spincolor*dirac
  inline void unsafe_spincolor_prod_dirac(spincolor out,spincolor in,dirac_matr *m)
  {spincolor_put_to_zero(out);for(size_t id1=0;id1<4;id1++) color_summ_the_prod_complex(out[m->pos[id1]],in[id1],m->entr[id1]);}
  
  //dirac*spincolor
  inline void safe_dirac_prod_spincolor(spincolor out,dirac_matr *m,spincolor in)
  {spincolor tmp;unsafe_dirac_prod_spincolor(tmp,m,in);spincolor_copy(out,tmp);}
  
  //spincolor*dirac
  inline void safe_spincolor_prod_dirac(spincolor out,spincolor in,dirac_matr *m)
  {spincolor tmp;unsafe_spincolor_prod_dirac(tmp,in,m);spincolor_copy(out,tmp);}
  
  //su3*spincolor
  inline void unsafe_su3_prod_spincolor(spincolor out,su3 U,spincolor in)
  {for(size_t is=0;is<4;is++) unsafe_su3_prod_color(out[is],U,in[is]);}
  inline void su3_summ_the_prod_spincolor(spincolor out,su3 U,spincolor in)
  {for(size_t is=0;is<4;is++) su3_summ_the_prod_color(out[is],U,in[is]);}
  inline void su3_subt_the_prod_spincolor(spincolor out,su3 U,spincolor in)
  {for(size_t is=0;is<4;is++) su3_subt_the_prod_color(out[is],U,in[is]);}
  
  //su3^*spincolor
  inline void unsafe_su3_dag_prod_spincolor(spincolor out,su3 U,spincolor in)
  {for(size_t is=0;is<4;is++) unsafe_su3_dag_prod_color(out[is],U,in[is]);}
  inline void safe_su3_dag_prod_spincolor(spincolor out,su3 U,spincolor in)
  {spincolor temp;unsafe_su3_dag_prod_spincolor(temp,U,in);spincolor_copy(out,temp);}
  
  inline void su3_dag_summ_the_prod_spincolor(spincolor out,su3 U,spincolor in)
  {for(size_t is=0;is<4;is++) su3_dag_summ_the_prod_color(out[is],U,in[is]);}
  inline void su3_dag_subt_the_prod_spincolor(spincolor out,su3 U,spincolor in)
  {for(size_t is=0;is<4;is++) su3_dag_subt_the_prod_color(out[is],U,in[is]);}
  
  //su3^*gamma*spincolor
  inline void unsafe_su3_dag_dirac_prod_spincolor(spincolor out,su3 U,dirac_matr *m,spincolor in)
  {
    color tmp;
    for(size_t id1=0;id1<4;id1++)
      {
	for(size_t ic=0;ic<NCOL;ic++) unsafe_complex_prod(tmp[ic],m->entr[id1],in[m->pos[id1]][ic]);
	unsafe_su3_dag_prod_color(out[id1],U,tmp);
      }
  }
  
  inline void unsafe_su3_dag_dirac_summ_the_prod_spincolor(spincolor out,su3 U,dirac_matr *m,spincolor in)
  {
    color tmp;
    for(size_t id1=0;id1<4;id1++)
      {
	for(size_t ic=0;ic<NCOL;ic++) unsafe_complex_prod(tmp[ic],m->entr[id1],in[m->pos[id1]][ic]);
	su3_dag_summ_the_prod_color(out[id1],U,tmp);
      }
  }
  
  //su3*dirac*spincolor
  inline void unsafe_su3_dirac_prod_spincolor(spincolor out,su3 U,dirac_matr *m,spincolor in)
  {
    color tmp;
    for(size_t id1=0;id1<4;id1++)
      {
	for(size_t ic=0;ic<NCOL;ic++) unsafe_complex_prod(tmp[ic],m->entr[id1],in[m->pos[id1]][ic]);
	unsafe_su3_prod_color(out[id1],U,tmp);
      }
  }
  
  inline void unsafe_su3_dirac_subt_the_prod_spincolor(spincolor out,su3 U,dirac_matr *m,spincolor in)
  {
    color tmp;
    for(size_t id1=0;id1<4;id1++)
      {
	for(size_t ic=0;ic<NCOL;ic++) unsafe_complex_prod(tmp[ic],m->entr[id1],in[m->pos[id1]][ic]);
	su3_subt_the_prod_color(out[id1],U,tmp);
      }
  }
  
  ///////////////////////////////// su3*colorspinspin ///////////////////////////////////
  
  inline void unsafe_su3_prod_colorspinspin(colorspinspin a,su3 b,colorspinspin c)
  {
    for(size_t id_so=0;id_so<4;id_so++)
      for(size_t id_si=0;id_si<4;id_si++)
	for(size_t c1=0;c1<NCOL;c1++)
	  {
	    unsafe_complex_prod(a[c1][id_si][id_so],b[c1][0],c[0][id_si][id_so]);
	    for(size_t c2=1;c2<NCOL;c2++) complex_summ_the_prod(a[c1][id_si][id_so],b[c1][c2],c[c2][id_si][id_so]);
	  }
  }
  
  inline void unsafe_su3_dag_prod_colorspinspin(colorspinspin a,su3 b,colorspinspin c)
  {
    for(size_t id_so=0;id_so<4;id_so++)
      for(size_t id_si=0;id_si<4;id_si++)
	for(size_t c1=0;c1<NCOL;c1++)
	  {
	    unsafe_complex_conj1_prod(a[c1][id_si][id_so],b[0][c1],c[0][id_si][id_so]);
	    for(size_t c2=1;c2<NCOL;c2++)
	      complex_summ_the_conj1_prod(a[c1][id_si][id_so],b[c2][c1],c[c2][id_si][id_so]);
	  }
  }
  
  inline void su3_summ_the_prod_colorspinspin(colorspinspin a,su3 b,colorspinspin c)
  {
    for(size_t id_so=0;id_so<4;id_so++)
      for(size_t id_si=0;id_si<4;id_si++)
	for(size_t c1=0;c1<NCOL;c1++)
	  for(size_t c2=0;c2<NCOL;c2++)
	    complex_summ_the_prod(a[c1][id_si][id_so],b[c1][c2],c[c2][id_si][id_so]);
  }
  
  inline void su3_dag_subt_the_prod_colorspinspin(colorspinspin a,su3 b,colorspinspin c)
  {
    for(size_t id_so=0;id_so<4;id_so++)
      for(size_t id_si=0;id_si<4;id_si++)
	for(size_t c1=0;c1<NCOL;c1++)
	  for(size_t c2=0;c2<NCOL;c2++)
	    complex_subt_the_conj1_prod(a[c1][id_si][id_so],b[c2][c1],c[c2][id_si][id_so]);
  }
  
  inline void su3_dag_summ_the_prod_colorspinspin(colorspinspin a,su3 b,colorspinspin c)
  {
    for(size_t id_so=0;id_so<4;id_so++)
      for(size_t id_si=0;id_si<4;id_si++)
	for(size_t c1=0;c1<NCOL;c1++)
	  for(size_t c2=0;c2<NCOL;c2++)
	    complex_summ_the_conj1_prod(a[c1][id_si][id_so],b[c2][c1],c[c2][id_si][id_so]);
  }
  
  ///////////////////////////////////// colorspinspin ////////////////////////////////////
  
  //colorspinspin*real or complex
  inline void colorspinspin_prod_double(colorspinspin out,colorspinspin in,double factor)
  {for(size_t i=0;i<NCOL;i++) spinspin_prod_double(out[i],in[i],factor);}
  inline void colorspinspin_prodassign_double(colorspinspin c,double f) {colorspinspin_prod_double(c,c,f);}
  inline void colorspinspin_prodassign_idouble(colorspinspin c,double f)
  {for(uint32_t ic=0;ic<NCOL;ic++) spinspin_prodassign_idouble(c[ic],f);}
  inline void colorspinspin_prod_complex(colorspinspin out,colorspinspin in,complex factor)
  {for(size_t ic=0;ic<NCOL;ic++) for(size_t id=0;id<4;id++)for(size_t jd=0;jd<4;jd++)safe_complex_prod(out[ic][id][jd],in[ic][id][jd],factor);}
  inline void colorspinspin_prod_complex_conj(colorspinspin out,colorspinspin in,complex factor)
  {complex temp;complex_conj(temp,factor);colorspinspin_prod_complex(out,in,temp);}
  inline void colorspinspin_prodassign_complex(colorspinspin c,complex f)
  {colorspinspin_prod_complex(c,c,f);}
  inline void colorspinspin_prodassign_complex_conj(colorspinspin c,complex f)
  {colorspinspin_prod_complex_conj(c,c,f);}
  
  //colorspinspin summ
  inline void colorspinspin_summ(colorspinspin out,colorspinspin in1,colorspinspin in2) {for(size_t i=0;i<NCOL;i++) spinspin_summ(out[i],in1[i],in2[i]);}
  inline void colorspinspin_summassign(colorspinspin out,colorspinspin in) {colorspinspin_summ(out,out,in);}
  
  //colorspinspin subt
  inline void colorspinspin_subt(colorspinspin out,colorspinspin in1,colorspinspin in2) {for(size_t i=0;i<NCOL;i++) spinspin_subt(out[i],in1[i],in2[i]);}
  inline void colorspinspin_subtassign(colorspinspin out,colorspinspin in) {colorspinspin_subt(out,out,in);}
  
  //summ two colorspinspin with a factor
  inline void colorspinspin_summ_the_prod_double(colorspinspin out,colorspinspin in,double factor)
  {for(size_t i=0;i<NCOL;i++) spinspin_summ_the_prod_double(out[i],in[i],factor);}
  inline void colorspinspin_summ_the_prod_idouble(colorspinspin out,colorspinspin in,double factor)
  {for(size_t i=0;i<NCOL;i++) spinspin_summ_the_prod_idouble(out[i],in[i],factor);}
  
  //colorspinspin*complex
  inline void unsafe_colorspinspin_prod_complex(colorspinspin out,colorspinspin in,complex factor)
  {for(size_t i=0;i<NCOL;i++) unsafe_spinspin_prod_complex(out[i],in[i],factor);}
  
  /////////////////////////////// su3spinspin /////////////////////////////////////////////
  
  //colorspinspin*real or complex
  inline void su3spinspin_prod_double(su3spinspin out,su3spinspin in,double factor)
  {for(size_t i=0;i<NCOL;i++) colorspinspin_summ_the_prod_double(out[i],in[i],factor);}
  inline void su3spinspin_prodassign_double(su3spinspin out,double factor) {for(size_t i=0;i<NCOL;i++) colorspinspin_prodassign_double(out[i],factor);}
  inline void su3spinspin_prodassign_idouble(su3spinspin out,double factor) {for(size_t i=0;i<NCOL;i++) colorspinspin_prodassign_idouble(out[i],factor);}
  inline void su3spinspin_prod_complex(su3spinspin out,su3spinspin in,complex factor) {for(size_t i=0;i<NCOL;i++) colorspinspin_prod_complex(out[i],in[i],factor);}
  inline void su3spinspin_prod_complex_conj(su3spinspin out,su3spinspin in,complex factor)
  {for(size_t i=0;i<NCOL;i++) colorspinspin_prod_complex_conj(out[i],in[i],factor);}
  inline void su3spinspin_prodassign_complex(su3spinspin c,complex f){su3spinspin_prod_complex(c,c,f);}
  inline void su3spinspin_prodassign_complex_conj(su3spinspin c,complex f) {su3spinspin_prod_complex_conj(c,c,f);}
  
  //su3spinspin summ
  inline void su3spinspin_summ(su3spinspin out,su3spinspin in1,su3spinspin in2) {for(size_t i=0;i<NCOL;i++) colorspinspin_summ(out[i],in1[i],in2[i]);}
  inline void su3spinspin_summassign(su3spinspin out,su3spinspin in) {su3spinspin_summ(out,out,in);}
  
  //su3spinspin subt
  inline void su3spinspin_subt(su3spinspin out,su3spinspin in1,su3spinspin in2) {for(size_t i=0;i<NCOL;i++) colorspinspin_subt(out[i],in1[i],in2[i]);}
  inline void su3spinspin_subtassign(su3spinspin out,su3spinspin in) {su3spinspin_subt(out,out,in);}
  
  //summ two su3spinspin with a factor
  inline void su3spinspin_summ_the_prod_double(su3spinspin out,su3spinspin in,double factor)
  {for(size_t i=0;i<NCOL;i++) colorspinspin_summ_the_prod_double(out[i],in[i],factor);}
  inline void su3spinspin_summ_the_prod_idouble(su3spinspin out,su3spinspin in,double factor)
  {for(size_t i=0;i<NCOL;i++) colorspinspin_summ_the_prod_idouble(out[i],in[i],factor);}
  
  ////////////////////////////////////// su3spinspin /////////////////////////////////////
  
  //su3spinspin*complex
  inline void unsafe_su3spinspin_prod_complex(su3spinspin out,su3spinspin in,complex factor)
  {for(size_t i=0;i<NCOL;i++) unsafe_colorspinspin_prod_complex(out[i],in[i],factor);}
  inline void safe_su3spinspin_prod_complex(su3spinspin out,su3spinspin in,complex factor)
  {su3spinspin temp;su3spinspin_copy(temp,in);unsafe_su3spinspin_prod_complex(out,temp,factor);}
  
  inline void unsafe_su3_prod_su3spinspin(su3spinspin a,su3 b,su3spinspin c)
  {
    for(size_t id_so=0;id_so<4;id_so++)
      for(size_t id_si=0;id_si<4;id_si++)
	for(size_t c1=0;c1<NCOL;c1++)
	  for(size_t c3=0;c3<NCOL;c3++)
	    {
	      unsafe_complex_prod(a[c1][c3][id_si][id_so],b[c1][0],c[0][c3][id_si][id_so]);
	      for(size_t c2=1;c2<NCOL;c2++) complex_summ_the_prod(a[c1][c3][id_si][id_so],b[c1][c2],c[c2][c3][id_si][id_so]);
	    }
  }
  
  inline void unsafe_su3_dag_prod_su3spinspin(su3spinspin a,su3 b,su3spinspin c)
  {
    for(size_t id_so=0;id_so<4;id_so++)
      for(size_t id_si=0;id_si<4;id_si++)
	for(size_t c1=0;c1<NCOL;c1++)
	  for(size_t c3=0;c3<NCOL;c3++)
	    {
	      unsafe_complex_conj1_prod(a[c1][c3][id_si][id_so],b[0][c1],c[0][c3][id_si][id_so]);
	      for(size_t c2=1;c2<NCOL;c2++)
		complex_summ_the_conj1_prod(a[c1][c3][id_si][id_so],b[c2][c1],c[c2][c3][id_si][id_so]);
	    }
  }
  
  inline void su3_summ_the_prod_su3spinspin(su3spinspin a,su3 b,su3spinspin c)
  {
    for(size_t id_so=0;id_so<4;id_so++)
      for(size_t id_si=0;id_si<4;id_si++)
	for(size_t c1=0;c1<NCOL;c1++)
	  for(size_t c3=0;c3<NCOL;c3++)
	    for(size_t c2=0;c2<NCOL;c2++)
	      complex_summ_the_prod(a[c1][c3][id_si][id_so],b[c1][c2],c[c2][c3][id_si][id_so]);
  }
  
  inline void su3_dag_subt_the_prod_su3spinspin(su3spinspin a,su3 b,su3spinspin c)
  {
    for(size_t id_so=0;id_so<4;id_so++)
      for(size_t id_si=0;id_si<4;id_si++)
	for(size_t c1=0;c1<NCOL;c1++)
	  for(size_t c2=0;c2<NCOL;c2++)
	    for(size_t c3=0;c3<NCOL;c3++)
	      complex_subt_the_conj1_prod(a[c1][c3][id_si][id_so],b[c2][c1],c[c2][c3][id_si][id_so]);
  }
  
  inline void su3_dag_summ_the_prod_su3spinspin(su3spinspin a,su3 b,su3spinspin c)
  {
    for(size_t id_so=0;id_so<4;id_so++)
      for(size_t id_si=0;id_si<4;id_si++)
	for(size_t c1=0;c1<NCOL;c1++)
	  for(size_t c2=0;c2<NCOL;c2++)
	    for(size_t c3=0;c3<NCOL;c3++)
	      complex_summ_the_conj1_prod(a[c1][c3][id_si][id_so],b[c2][c1],c[c2][c3][id_si][id_so]);
  }
  
  //put a matrix to random used passed random generator
  inline void su3_put_to_rnd(su3 u_ran,rnd_gen &rnd)
  {
    su3_put_to_id(u_ran);
    
    for(size_t i1=0;i1<NCOL;i1++)
      for(size_t i2=i1+1;i2<NCOL;i2++)
	{
	  //generate u0,u1,u2,u3 random on the four dim. sphere
	  double u0=rnd_get_unif(&rnd,-1,1);
	  double alpha=sqrt(1-u0*u0);
	  double phi=rnd_get_unif(&rnd,0,2*M_PI);
	  double costheta=rnd_get_unif(&rnd,-1,1);
	  double sintheta=sqrt(1-costheta*costheta);
	  double u3=alpha*costheta;
	  double u1=alpha*sintheta*cos(phi);
	  double u2=alpha*sintheta*sin(phi);
	  
	  //define u_l as unit matrix ...
	  su3 u_l;
	  su3_put_to_id(u_l);
	  
	  //... and then modify the elements in the chosen su(2) subgroup
	  u_l[i1][i1][RE]=u0;
	  u_l[i1][i1][IM]=u3;
	  u_l[i1][i2][RE]=u2;
	  u_l[i1][i2][IM]=u1;
	  u_l[i2][i1][RE]=-u2;
	  u_l[i2][i1][IM]=u1;
	  u_l[i2][i2][RE]=u0;
	  u_l[i2][i2][IM]=-u3;
	  
	  safe_su3_prod_su3(u_ran,u_l,u_ran);
	}
  }
  
  void anti_hermitian_exact_i_exponentiate_ingredients(anti_hermitian_exp_ingredients &out,su3 Q);
  
  //build the exponential from the ingredients
  inline void safe_anti_hermitian_exact_i_exponentiate(su3 out,anti_hermitian_exp_ingredients &ing)
  {
    CRASH_IF_NOT_3COL();
    
    //compute out according to (eq. 13)
    su3_put_to_diag(out,ing.f[0]);
    su3_summ_the_prod_complex(out,ing.Q,ing.f[1]);
    su3_summ_the_prod_complex(out,ing.Q2,ing.f[2]);
  }
  inline void safe_anti_hermitian_exact_i_exponentiate(su3 out,su3 Q)
  {
    anti_hermitian_exp_ingredients ing;
    anti_hermitian_exact_i_exponentiate_ingredients(ing,Q);
    safe_anti_hermitian_exact_i_exponentiate(out,ing);
  }
  
  //can be used for an hermitian matrix
  inline void safe_hermitian_exact_exponentiate(su3 out,su3 in)
  {
    su3 Q;
    su3_prod_idouble(Q,in,-1);
    
    safe_anti_hermitian_exact_i_exponentiate(out,Q);
  }
  
  //return sqrt(|U*U^+-1|)
  inline double su3_get_non_unitariness(su3 u)
  {
    su3 zero;
    su3_put_to_id(zero);
    su3_subt_the_prod_su3_dag(zero,u,u);
    
    return sqrt(real_part_of_trace_su3_prod_su3_dag(zero,zero));
  }
}

#endif
