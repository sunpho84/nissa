#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <string.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/random.hpp"

#include "complex.hpp"
#include "float_128.hpp"
#include "new_types_definitions.hpp"
#include "su3.hpp"

#include "operations/su3_paths/squared_staples.hpp"
#include "routines/ios.hpp"
#include "routines/math_routines.hpp"

namespace nissa
{
  void su3_put_to_zero(su3 m) {memset(m,0,sizeof(su3));}
  void as2t_su3_put_to_zero(as2t_su3 m) {memset(m,0,sizeof(as2t_su3));}
  void spincolor_put_to_zero(spincolor m) {memset(m,0,sizeof(spincolor));}
  void su3spinspin_put_to_zero(su3spinspin m) {memset(m,0,sizeof(su3spinspin));}
  void su3_put_to_id(su3 m) {su3_put_to_zero(m);for(int ic=0;ic<3;ic++) m[ic][ic][0]=1;}
  void su3_put_to_diag(su3 m,color in) {su3_put_to_zero(m);for(int ic=0;ic<3;ic++) complex_copy(m[ic][ic],in[ic]);}
  void su3_put_to_diag(su3 m,complex in) {su3_put_to_zero(m);for(int ic=0;ic<3;ic++) complex_copy(m[ic][ic],in);}
  void su3_put_to_diag(su3 m,double in) {su3_put_to_zero(m);for(int ic=0;ic<3;ic++) m[ic][ic][0]=in;}
  
  //////////////////////////////////////// Copy /////////////////////////////////////
  
  void color_copy(color b,color a) {memcpy(b,a,sizeof(color));}
  void su3_copy(su3 b,su3 a) {memcpy(b,a,sizeof(su3));}
  void quad_su3_copy(quad_su3 b,quad_su3 a) {memcpy(b,a,sizeof(quad_su3));}
  void spincolor_copy(spincolor b,spincolor a) {memcpy(b,a,sizeof(spincolor));}
  void colorspinspin_copy(colorspinspin b,colorspinspin a) {memcpy(b,a,sizeof(colorspinspin));}
  
  //////////////////// Switch directions so to agree to ILDG ordering ////////////////////
  
  void quad_su3_nissa_to_ildg_reord(quad_su3 out,quad_su3 in)
  {
    quad_su3 buff;
    quad_su3_copy(buff,in);
    
    su3_copy(out[3],buff[0]);
    su3_copy(out[0],buff[1]);
    su3_copy(out[1],buff[2]);
    su3_copy(out[2],buff[3]);
  }
  
  void quad_su3_ildg_to_nissa_reord(quad_su3 out,quad_su3 in)
  {
    quad_su3 buff;
    quad_su3_copy(buff,in);
    
    su3_copy(out[0],buff[3]);
    su3_copy(out[1],buff[0]);
    su3_copy(out[2],buff[1]);
    su3_copy(out[3],buff[2]);
  }
  
  ////////////////////////////////// Operations between colors //////////////////////////
  
  //just print a color
  void color_print(color c)
  {
    for(int ic=0;ic<3;ic++) printf("%+016.16le,%+016.16le\t",c[ic][0],c[ic][1]);
    master_printf("\n");
  }
  
  //summ two colors
  void color_summ(color a,color b,color c) {for(int i=0;i<6;i++) ((double*)a)[i]=((double*)b)[i]+((double*)c)[i];}
  
  void color_isumm(color a,color b,color c)
  {for(int i=0;i<6;i+=2) {((double*)a)[i]=((double*)b)[i]-((double*)c)[i+1];((double*)a)[i+1]=((double*)b)[i+1]+((double*)c)[i];}}
  
  void color_isubt(color a,color b,color c)
  {for(int i=0;i<6;i+=2) {((double*)a)[i]=((double*)b)[i]+((double*)c)[i+1];((double*)a)[i+1]=((double*)b)[i+1]-((double*)c)[i];}}
  
  void color_subt(color a,color b,color c) {for(int i=0;i<6;i++) ((double*)a)[i]=((double*)b)[i]-((double*)c)[i];}
  void color_summassign(color a,color b) {for(int i=0;i<6;i++) ((double*)a)[i]+=((double*)b)[i];}
  void color_subtassign(color a,color b) {for(int i=0;i<6;i++) ((double*)a)[i]-=((double*)b)[i];}
  
  void color_isummassign(color a,color b) {for(int i=0;i<6;i+=2) {((double*)a)[i]-=((double*)b)[i+1];((double*)a)[i+1]+=((double*)b)[i];}}
  
  void color_isubtassign(color a,color b)
  {for(int i=0;i<6;i+=2) {((double*)a)[i]+=((double*)b)[i+1];((double*)a)[i+1]-=((double*)b)[i];}}
  
  void color_prod_double(color a,color b,double c)
  {for(int i=0;i<6;i++) ((double*)a)[i]=((double*)b)[i]*c;}
  
  void color_scalar_prod(complex out,color a,color b)
  {
    unsafe_complex_conj2_prod(out,a[0],b[0]);
    complex_summ_the_conj2_prod(out,a[1],b[1]);
    complex_summ_the_conj2_prod(out,a[2],b[2]);
  }
  
  /////////////////////////////// Generate an hermitean matrix ///////////////////////
  
  //Taken from M.D'Elia
  void herm_put_to_gauss(su3 H,rnd_gen *gen,double sigma)
  {
    const double one_by_sqrt3=0.577350269189626;
    const double two_by_sqrt3=1.15470053837925;
    
    double r[8];
    for(int ir=0;ir<4;ir++)
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
  
  // A gauss vector has complex components z which are gaussian distributed
  // with <z~ z> = sigma
  void color_put_to_gauss(color H,rnd_gen *gen,double sigma)
  {
    complex ave={0,0};
    rnd_get_gauss_complex(H[0],gen,ave,sigma);
    rnd_get_gauss_complex(H[1],gen,ave,sigma);
    rnd_get_gauss_complex(H[2],gen,ave,sigma);
  }
  
  ////////////////////////////////// Operations between su3 //////////////////////////
  
  //comput third row according to other 2
  void su3_build_third_row(su3 o)
  {
    unsafe_complex_conj_conj_prod(o[2][0],o[0][1],o[1][2]);
    complex_subt_the_conj_conj_prod(o[2][0],o[0][2],o[1][1]);
    unsafe_complex_conj_conj_prod(o[2][1],o[0][2],o[1][0]);
    complex_subt_the_conj_conj_prod(o[2][1],o[0][0],o[1][2]);
    unsafe_complex_conj_conj_prod(o[2][2],o[0][0],o[1][1]);
    complex_subt_the_conj_conj_prod(o[2][2],o[0][1],o[1][0]);
  }
  
  //just print an su3 matrix
  void su3_print(su3 U)
  {
    for(int ic1=0;ic1<3;ic1++)
      {
	for(int ic2=0;ic2<3;ic2++) printf("%+016.16le,%+016.16le\t",U[ic1][ic2][0],U[ic1][ic2][1]);
	printf("\n");
      }
    printf("\n");
  }
  
  //return the trace of an su3 matrix
  void su3_trace(complex tr,su3 m)
  {
    tr[0]=m[0][0][0];
    tr[1]=m[0][0][1];
    
    complex_summ(tr,tr,m[1][1]);
    complex_summ(tr,tr,m[2][2]);
  }
  
  //return only the real part of an su3 matrix
  double su3_real_trace(su3 m)
  {return m[0][0][RE]+m[1][1][RE]+m[2][2][RE];}
  
  //take projection of the su2 matrix over an su3 matrix
  //return the inverse modulo of the part parallel in the original matrix
  double su2_part_of_su3(double &A,double &B,double &C,double &D,su3 in,int isub_gr)
  {
    //take indices of the subgroup
    int a=su3_sub_gr_indices[isub_gr][0];
    int b=su3_sub_gr_indices[isub_gr][1];
    
    //extract part parallel to sigmas
    A=in[a][a][0]+in[b][b][0];
    B=in[a][b][1]+in[b][a][1];
    C=in[a][b][0]-in[b][a][0];
    D=in[a][a][1]-in[b][b][1];
    
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
  double su2_part_of_su3(su2 out,su3 in,int isub_gr)
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
  void su2_prodassign_su3(su2 mod,int isub_gr,su3 in)
  {
    int ic1=su3_sub_gr_indices[isub_gr][0];
    int ic2=su3_sub_gr_indices[isub_gr][1];
    
    //create the two new rows of the matrix, column by column
    for(int ic=0;ic<3;ic++)
      {
	//first row
	complex row1;
	safe_complex_prod    (row1,mod[0][0],in[ic1][ic]);
	complex_summ_the_prod(row1,mod[0][1],in[ic2][ic]);
	//second row
	complex row2;      
	safe_complex_prod    (row2,mod[1][0],in[ic1][ic]);
	complex_summ_the_prod(row2,mod[1][1],in[ic2][ic]);
	
	//change the two lines in the matrix
	complex_copy(in[ic1][ic],row1);
	complex_copy(in[ic2][ic],row2);
      }
  }
  
  //in the form A+B*i*sigma1+...
  void su2_prodassign_su3(double A,double B,double C,double D,int isub_gr,su3 in)
  {
    su2 mod={{{A,D},{C,B}},{{-C,B},{A,-D}}};
    su2_prodassign_su3(mod,isub_gr,in);
  }
  
  //summ the trace to the input
  void su3_summ_the_trace(complex tr,su3 m)
  {
    complex_summ(tr,tr,m[0][0]);
    complex_summ(tr,tr,m[1][1]);
    complex_summ(tr,tr,m[2][2]);
  }
  
  //return the anti-hermitian traceless part of an su3 matrix
  void unsafe_su3_traceless_anti_hermitian_part(su3 out,su3 in)
  {
    double trace_im_third=(in[0][0][1]+in[1][1][1]+in[2][2][1])/3;
    
    //real part of diagonal: 0
    out[0][0][0]=out[1][1][0]=out[2][2][0]=0;
    //imag part of diagonal: subtract the trace
    out[0][0][1]=in[0][0][1]-trace_im_third;
    out[1][1][1]=in[1][1][1]-trace_im_third;
    out[2][2][1]=in[2][2][1]-trace_im_third;
    //out-of-diag real part
    out[1][0][0]=-(out[0][1][0]=(in[0][1][0]-in[1][0][0])/2);
    out[2][0][0]=-(out[0][2][0]=(in[0][2][0]-in[2][0][0])/2);
    out[2][1][0]=-(out[1][2][0]=(in[1][2][0]-in[2][1][0])/2);
    //out-of-diag imag part
    out[1][0][1]=out[0][1][1]=(in[0][1][1]+in[1][0][1])/2;
    out[2][0][1]=out[0][2][1]=(in[0][2][1]+in[2][0][1])/2;
    out[2][1][1]=out[1][2][1]=(in[1][2][1]+in[2][1][1])/2;
  }
  
  //return the hermitian traceless part of an su3 matrix
  void unsafe_su3_traceless_hermitian_part(su3 out,su3 in)
  {
    double trace_re_third=(in[0][0][0]+in[1][1][0]+in[2][2][0])/3;
    
    //imag part of diagonal: 0
    out[0][0][1]=out[1][1][1]=out[2][2][1]=0;
    //real part of diagonal: subtract the trace
    out[0][0][0]=in[0][0][0]-trace_re_third;
    out[1][1][0]=in[1][1][0]-trace_re_third;
    out[2][2][0]=in[2][2][0]-trace_re_third;
    //out-of-diag real part
    out[1][0][0]=out[0][1][0]=(in[0][1][0]+in[1][0][0])/2;
    out[2][0][0]=out[0][2][0]=(in[0][2][0]+in[2][0][0])/2;
    out[2][1][0]=out[1][2][0]=(in[1][2][0]+in[2][1][0])/2;
    //out-of-diag imag part
    out[1][0][1]=-(out[0][1][1]=(in[0][1][1]-in[1][0][1])/2);
    out[2][0][1]=-(out[0][2][1]=(in[0][2][1]-in[2][0][1])/2);
    out[2][1][1]=-(out[1][2][1]=(in[1][2][1]-in[2][1][1])/2);
  }
  
  //calculate the determinant of an su3 matrix
  void su3_det(complex d,su3 U)
  {
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
  }
  
  //calculate the real part of the determinant of an su3 matrix
  double su3_real_det(su3 u)
  {
    return
      u[0][2][IM]*(u[1][1][RE]*u[2][0][IM]+u[1][1][IM]*u[2][0][RE]-u[1][0][RE]*u[2][1][IM]-u[1][0][IM]*u[2][1][RE])+
      u[0][2][RE]*(u[1][1][IM]*u[2][0][IM]-u[1][1][RE]*u[2][0][RE]-u[1][0][IM]*u[2][1][IM]+u[1][0][RE]*u[2][1][RE])+
      u[0][1][IM]*(u[1][0][RE]*u[2][2][IM]-u[1][2][RE]*u[2][0][IM]-u[1][2][IM]*u[2][0][RE]+u[1][0][IM]*u[2][2][RE])+
      u[0][1][RE]*(u[1][0][IM]*u[2][2][IM]-u[1][2][IM]*u[2][0][IM]+u[1][2][RE]*u[2][0][RE]-u[1][0][RE]*u[2][2][RE])+
      u[0][0][IM]*(u[1][2][RE]*u[2][1][IM]+u[1][2][IM]*u[2][1][RE]-u[1][1][RE]*u[2][2][IM]-u[1][1][IM]*u[2][2][RE])+
      u[0][0][RE]*(u[1][2][IM]*u[2][1][IM]-u[1][2][RE]*u[2][1][RE]-u[1][1][IM]*u[2][2][IM]+u[1][1][RE]*u[2][2][RE]);
  }
  
  //return the hemitian su3 matrix
  void unsafe_su3_hermitian(su3 out,su3 in)
  {
    for(int ic_in=0;ic_in<3;ic_in++)
      for(int ic_out=0;ic_out<3;ic_out++)
	{
	  out[ic_in][ic_out][0]= in[ic_out][ic_in][0];
	  out[ic_in][ic_out][1]=-in[ic_out][ic_in][1];
	}
  }
  void safe_su3_hermitian(su3 out,su3 in)
  {
    su3 tmp;
    unsafe_su3_hermitian(tmp,in);
    su3_copy(out,tmp);
  }
  
  //return the transposed su3 matrix
  void unsafe_su3_transpose(su3 out,su3 in)
  {
    for(int ic_in=0;ic_in<3;ic_in++)
      for(int ic_out=0;ic_out<3;ic_out++)
	{
	  out[ic_in][ic_out][0]=in[ic_out][ic_in][0];
	  out[ic_in][ic_out][1]=in[ic_out][ic_in][1];
	}
  }
  void safe_su3_transpose(su3 out,su3 in)
  {
    su3 tmp;
    unsafe_su3_transpose(tmp,in);
    su3_copy(out,tmp);
  }
  
  //summ two su3 matrixes
  void su3_summ(su3 a,su3 b,su3 c) {for(int i=0;i<18;i++) ((double*)a)[i]=((double*)b)[i]+((double*)c)[i];}
  void su3_summassign(su3 a,su3 b){su3_summ(a,a,b);}
  void su3_summ_real(su3 a,su3 b,double c) {su3_copy(a,b);for(int i=0;i<3;i++) a[i][i][0]=b[i][i][0]+c;}
  void su3_subt(su3 a,su3 b,su3 c) {for(int i=0;i<18;i++) ((double*)a)[i]=((double*)b)[i]-((double*)c)[i];}
  void su3_subt_complex(su3 a,su3 b,complex c)
  {for(int i=0;i<3;i++) for(int ri=0;ri<2;ri++) a[i][i][ri]=b[i][i][ri]-c[ri];}
  void unsafe_su3_subt_su3_dag(su3 a,su3 b,su3 c)
  {
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	complex_subt_conj2(a[i][j],b[i][j],c[j][i]);
  }
  
  //Product of two su3 matrixes
  void unsafe_su3_prod_su3(su3 a,su3 b,su3 c,int nr_max)
  {
    for(int ir_out=0;ir_out<nr_max;ir_out++)
      for(int ic_out=0;ic_out<3;ic_out++)
	{
	  unsafe_complex_prod(a[ir_out][ic_out],b[ir_out][0],c[0][ic_out]);
	  for(int itemp=1;itemp<3;itemp++)
	    complex_summ_the_prod(a[ir_out][ic_out],b[ir_out][itemp],c[itemp][ic_out]);
	}
  }
  void safe_su3_prod_su3(su3 a,su3 b,su3 c) {su3 d;unsafe_su3_prod_su3(d,b,c);su3_copy(a,d);}
  void su3_summ_the_prod_su3(su3 a,su3 b,su3 c)
  {
    for(int ir_out=0;ir_out<3;ir_out++)
      for(int ic_out=0;ic_out<3;ic_out++)
	for(int itemp=0;itemp<3;itemp++)
	  complex_summ_the_prod(a[ir_out][ic_out],b[ir_out][itemp],c[itemp][ic_out]);
  }
  void su3_summ_the_prod_su3_dag(su3 a,su3 b,su3 c)
  {
    for(int ir_out=0;ir_out<3;ir_out++)
      for(int ic_out=0;ic_out<3;ic_out++)
	for(int itemp=0;itemp<3;itemp++)
	  complex_summ_the_conj2_prod(a[ir_out][ic_out],b[ir_out][itemp],c[ic_out][itemp]);
  }
  void su3_summ_the_dag_prod_su3(su3 a,su3 b,su3 c)
  {
    for(int ir_out=0;ir_out<3;ir_out++)
      for(int ic_out=0;ic_out<3;ic_out++)
	for(int itemp=0;itemp<3;itemp++)
	  complex_summ_the_conj1_prod(a[ir_out][ic_out],b[itemp][ir_out],c[itemp][ic_out]);
  }
  
  //Product of two su3 matrixes
  void unsafe_su3_dag_prod_su3(su3 a,su3 b,su3 c,int nr_max)
  {
    for(int ir_out=0;ir_out<nr_max;ir_out++)
      for(int ic_out=0;ic_out<3;ic_out++)
	{
	  unsafe_complex_conj1_prod(a[ir_out][ic_out],b[0][ir_out],c[0][ic_out]);
	  for(int itemp=1;itemp<3;itemp++)
	    complex_summ_the_conj1_prod(a[ir_out][ic_out],b[itemp][ir_out],c[itemp][ic_out]);
	}
  }
  void safe_su3_dag_prod_su3(su3 a,su3 b,su3 c) {su3 d;unsafe_su3_dag_prod_su3(d,b,c);su3_copy(a,d);}
  void su3_dag_summ_the_prod_su3(su3 a,su3 b,su3 c)
  {
    for(int ir_out=0;ir_out<3;ir_out++)
      for(int ic_out=0;ic_out<3;ic_out++)
	for(int itemp=0;itemp<3;itemp++)
	  complex_summ_the_conj1_prod(a[ir_out][ic_out],b[itemp][ir_out],c[itemp][ic_out]);
  }
  
  //Product of two su3 matrixes
  void unsafe_su3_prod_su3_dag(su3 a,su3 b,su3 c,int nr_max)
  {
    for(int ir_out=0;ir_out<nr_max;ir_out++)
      for(int ic_out=0;ic_out<3;ic_out++)
	{
	  unsafe_complex_conj2_prod(a[ir_out][ic_out],b[ir_out][0],c[ic_out][0]);
	  complex_summ_the_conj2_prod(a[ir_out][ic_out],b[ir_out][1],c[ic_out][1]);
	  complex_summ_the_conj2_prod(a[ir_out][ic_out],b[ir_out][2],c[ic_out][2]);
	}
  }
  void safe_su3_prod_su3_dag(su3 a,su3 b,su3 c) {su3 d;unsafe_su3_prod_su3_dag(d,b,c);su3_copy(a,d);}
  
  //subtract the product
  void su3_subt_the_prod_su3_dag(su3 a,su3 b,su3 c)
  {
    for(int ir_out=0;ir_out<3;ir_out++)
      for(int ic_out=0;ic_out<3;ic_out++)
	{
	  complex_subt_the_conj2_prod(a[ir_out][ic_out],b[ir_out][0],c[ic_out][0]);
	  complex_subt_the_conj2_prod(a[ir_out][ic_out],b[ir_out][1],c[ic_out][1]);
	  complex_subt_the_conj2_prod(a[ir_out][ic_out],b[ir_out][2],c[ic_out][2]);
	}
  }
  
  //Trace of the product of two su3 matrices
  double real_part_of_trace_su3_prod_su3_dag(su3 a,su3 b)
  {
    double t=0;
    
    for(int ic1=0;ic1<3;ic1++)
      for(int ic2=0;ic2<3;ic2++)
	t+=a[ic1][ic2][0]*b[ic1][ic2][0]+a[ic1][ic2][1]*b[ic1][ic2][1];
    
    return t;
  }
  
  //Trace of the product of two su3 matrices
  void trace_su3_prod_su3(complex t,su3 a,su3 b)
  {
    t[RE]=t[IM]=0;
    for(int ic1=0;ic1<3;ic1++)
      for(int ic2=0;ic2<3;ic2++)
	complex_summ_the_prod(t,a[ic1][ic2],b[ic2][ic1]);
  }
  
  //Product of two su3 matrixes
  void unsafe_su3_dag_prod_su3_dag(su3 a,su3 b,su3 c,int nr_max)
  {
    for(int ir_out=0;ir_out<nr_max;ir_out++)
      for(int ic_out=0;ic_out<3;ic_out++)
	{
	  unsafe_complex_conj_conj_prod(a[ir_out][ic_out],b[0][ir_out],c[ic_out][0]);
	  for(int itemp=1;itemp<3;itemp++)
	    complex_summ_the_conj_conj_prod(a[ir_out][ic_out],b[itemp][ir_out],c[ic_out][itemp]);
	}
  }
  void safe_su3_dag_prod_su3_dag(su3 a,su3 b,su3 c)
  {su3 d;unsafe_su3_dag_prod_su3_dag(d,b,c);su3_copy(a,d);}
  
  //product of an su3 matrix by a complex
  void unsafe_su3_prod_complex(su3 a,su3 b,complex c)
  {
    complex *ca=(complex*)a;
    complex *cb=(complex*)b;
    
    for(int i=0;i<9;i++) unsafe_complex_prod(ca[i],cb[i],c);
  }
  void safe_su3_prod_complex(su3 a,su3 b,complex c)
  {su3 d;unsafe_su3_prod_complex(d,b,c);su3_copy(a,d);}
  
  //product of an su3 matrix by a complex
  void unsafe_su3_prod_conj_complex(su3 a,su3 b,complex c)
  {
    complex *ca=(complex*)a;
    complex *cb=(complex*)b;
    
    for(int i=0;i<9;i++) unsafe_complex_conj2_prod(ca[i],cb[i],c);
  }
  void safe_su3_prod_conj_complex(su3 a,su3 b,complex c)
  {su3 d;unsafe_su3_prod_conj_complex(d,b,c);su3_copy(a,d);}
  
  //product of an su3 matrix by a real
  void su3_prod_double(su3 a,su3 b,double r)
  {
    double *da=(double*)a;
    double *db=(double*)b;
    
    for(int i=0;i<18;i++) da[i]=db[i]*r;
  }
  
  //hermitian of su3 matrix times a real
  void unsafe_su3_hermitian_prod_double(su3 a,su3 b,double r)
  {
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	{
	  a[i][j][0]= r*b[j][i][0];
	  a[i][j][1]=-r*b[j][i][1];
	}
  }
  void safe_su3_hermitian_prod_double(su3 a,su3 b,double r)
  {
    for(int i=0;i<3;i++)
      {
	a[i][i][0]=+r*b[i][i][0];
	a[i][i][1]=-r*b[i][i][1];
	for(int j=i+1;j<3;j++)
	  {
	    double a_i_j_0=+r*b[j][i][0];
	    double a_i_j_1=-r*b[j][i][1];
	    a[j][i][0]=+r*b[i][j][0];
	    a[j][i][1]=-r*b[i][j][1];
	    a[i][j][0]=a_i_j_0;
	    a[i][j][1]=a_i_j_1;
	  }
      }
  }
  void su3_summ_the_hermitian_prod_double(su3 a,su3 b,double r)
  {
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	{
	  a[i][j][0]+= r*b[j][i][0];
	  a[i][j][1]+=-r*b[j][i][1];
	}
  }
  
  //summ the prod of su3 with imag
  void su3_prod_idouble(su3 a,su3 b,double r)
  {
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	{
	  double c=b[i][j][0]*r;
	  a[i][j][0]=-b[i][j][1]*r;
	  a[i][j][1]=c;
	}
  }
  void su3_summ_the_prod_idouble(su3 a,su3 b,double r)
  {
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	{
	  a[i][j][0]-=b[i][j][1]*r;
	  a[i][j][1]+=b[i][j][0]*r;
	}
  }
  
  //summ the prod of su3 with real
  void su3_summ_the_prod_double(su3 a,su3 b,double r)
  {
    double *da=(double*)a;
    double *db=(double*)b;
    
    for(int i=0;i<18;i++) da[i]+=db[i]*r;
  }
  
  //summ the prod of the dag su3 with real
  void su3_dag_summ_the_prod_double(su3 a,su3 b,double r)
  {
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	{
	  a[i][j][RE]+=b[j][i][RE]*r;
	  a[i][j][IM]-=b[j][i][IM]*r;
	}
  }
  
  //combine linearly two su3 elements
  void su3_linear_comb(su3 a,su3 b,double cb,su3 c,double cc)
  {for(int i=0;i<18;i++) ((double*)a)[i]=((double*)b)[i]*cb+((double*)c)[i]*cc;}
  
  //summ the prod of su3 with complex
  void su3_summ_the_prod_complex(su3 a,su3 b,complex c)
  {
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	complex_summ_the_prod(a[i][j],b[i][j],c);
  }
  
  //calculate explicitely the inverse
  void unsafe_su3_explicit_inverse(su3 invU,su3 U)
  {
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
    
    for(int icol1=0;icol1<3;icol1++)
      for(int icol2=0;icol2<3;icol2++)
	safe_complex_prod(invU[icol1][icol2],invU[icol1][icol2],rec_det);
  }
  void safe_su3_explicit_inverse(su3 invU,su3 U)
  {su3 tempU;unsafe_su3_explicit_inverse(tempU,U);su3_copy(invU,tempU);}
  
  //summ of the squared norm of the entries
  double su3_normq(su3 U)
  {
    double normq=0;
    
    for(int icol1=0;icol1<3;icol1++)
      for(int icol2=0;icol2<3;icol2++)
	for(int ri=0;ri<2;ri++)
	  normq+=U[icol1][icol2][ri]*U[icol1][icol2][ri];
    
    return normq;
  }
  
  //compute the square root of y numerically
  void unsafe_su3_sqrt(su3 x,su3 y)
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
	err=sqrt(su3_normq(t));
      }
    while(err>1.e-15);
  }
  void safe_su3_sqrt(su3 x,su3 y) {su3 t;unsafe_su3_sqrt(t,y);su3_copy(x,t);}
  
  //exponentiate an su3 matrix through taylor expansion
  void unsafe_su3_taylor_exponentiate(su3 out,su3 in,int order)
  {
    //1st terms
    double coef=1;
    su3 temp;
    su3_copy(temp,in);
    
    //order 0+1
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	{
	  out[i][j][0]=(i==j) ? in[i][j][0]+1 : in[i][j][0];
	  out[i][j][1]=in[i][j][1];
	}
    
    //higher orders
    for(int iorder=2;iorder<=order;iorder++)
      {
	safe_su3_prod_su3(temp,temp,in);
	coef/=iorder;
	su3_summ_the_prod_double(out,temp,coef);
      }
  }
  
  //unitarize by orthonormalizing the rows
  void su3_unitarize_orthonormalizing(su3 o,su3 i)
  {
    //compute the squared norm of row 0
    double row0_sq_norm=squared_complex_norm(i[0][0])+squared_complex_norm(i[0][1])+squared_complex_norm(i[0][2]);
    
    //compute the scalar product between row 1 and 0
    complex row10_sc_prod;
    unsafe_complex_conj2_prod(  row10_sc_prod,i[1][0],i[0][0]);
    complex_summ_the_conj2_prod(row10_sc_prod,i[1][1],i[0][1]);
    complex_summ_the_conj2_prod(row10_sc_prod,i[1][2],i[0][2]);
    
    //orthonormalize row 1 
    complex f;
    complex_prod_double(f,row10_sc_prod,1/row0_sq_norm);
    
    for(int c=0;c<3;c++)
      {
	for(int ri=0;ri<2;ri++)
	  o[1][c][ri]=i[1][c][ri];
	complex_subt_the_prod(o[1][c],f,i[0][c]);
      }
    
    double row0_norm=1/sqrt(row0_sq_norm);
    double row1_norm=1/sqrt(squared_complex_norm(o[1][0])+squared_complex_norm(o[1][1])+squared_complex_norm(o[1][2]));
    
    //normalize the rows
    for(int c=0;c<3;c++)
      for(int ri=0;ri<2;ri++)
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
  
  //overrelax the link
  void su3_overrelax(su3 out,su3 in,double w,double *coeff,const int ord)
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
  
  //overrelax the link using approximated exponentiation
  void su3_overrelax(su3 out,su3 in,double w)
  {
    double coeff[5]={1,w,w*(w-1)/2,w*(w-1)*(w-2)/6,w*(w-1)*(w-2)*(w-3)/24};
    su3_overrelax(out,in,w,coeff,5);
  }
  
  //unitarize an su3 matrix by taking explicitely the inverse and averaging with it
  void su3_unitarize_explicitly_inverting(su3 new_link,su3 prop_link)
  {
    su3 inv,temp_link;
    double gamma,residue;
    
    memcpy(temp_link,prop_link,sizeof(su3));
    
    do
      {
	unsafe_su3_explicit_inverse(inv,temp_link);
	gamma=sqrt(su3_normq(inv)/su3_normq(temp_link));
	
	//average U and U^-1^+
	residue=0;
	for(int icol1=0;icol1<3;icol1++)
	  for(int icol2=0;icol2<3;icol2++)
	    {
	      new_link[icol1][icol2][0]=0.5*(temp_link[icol1][icol2][0]*gamma+inv[icol2][icol1][0]/gamma);
	      new_link[icol1][icol2][1]=0.5*(temp_link[icol1][icol2][1]*gamma-inv[icol2][icol1][1]/gamma);
	      for(int ri=0;ri<2;ri++)
		{
		  double diff=new_link[icol1][icol2][ri]-temp_link[icol1][icol2][ri];
		  residue+=diff*diff;
		}
	    }
	
	memcpy(temp_link,new_link,sizeof(su3));
	
	residue=sqrt(residue); 
      }
    while(residue>1.e-15);
    
    //divide by third root of det
    complex det,fact;
    su3_det(det,new_link);
    complex_pow(fact,det,-1.0/3);
    safe_su3_prod_complex(new_link,new_link,fact);
  }
  
  //perform an iteration of maximal projection trace
  //by maximing the trace over three different subgroups of su3
  void su3_unitarize_maximal_trace_projecting_iteration(su3 U,su3 M)
  {
    //loop over the three subgroups
    for(int isub_gr=0;isub_gr<3;isub_gr++)
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
  
  //perform maximal projection trace up to reaching the machine precision
  void su3_unitarize_maximal_trace_projecting(su3 U,su3 M)
  {
    //initialize the guess
    su3_unitarize_explicitly_inverting(U,M);
    
    //compute initial trace
    double new_trace=real_part_of_trace_su3_prod_su3_dag(U,M);
    double old_trace,residue;
    
    //loop up to reach machine precision
    do
      {
	old_trace=new_trace;
	su3_unitarize_maximal_trace_projecting_iteration(U,M);
	new_trace=real_part_of_trace_su3_prod_su3_dag(U,M);
	residue=2*fabs(new_trace-old_trace)/(new_trace+old_trace);
      }
    while(residue>1.e-15);
  }
  
  //return a single link after the heatbath procedure
  //routines to be shrinked!
  void su3_find_heatbath(su3 out,su3 in,su3 staple,double beta,int nhb_hits,rnd_gen *gen)
  {
    //compute the original contribution to the action due to the given link 
    su3 prod;
    unsafe_su3_prod_su3_dag(prod,in,staple);
    
    //copy in link to out
    su3_copy(out,in);
    
    //iterate over heatbath hits
    for(int ihit=0;ihit<nhb_hits;ihit++)
      //scan all the three possible subgroups
      for(int isub_gr=0;isub_gr<3;isub_gr++)
	{
	  //take the part of the su3 matrix
	  double r0,r1,r2,r3;
	  double smod=su2_part_of_su3(r0,r1,r2,r3,prod,isub_gr);
	  
	  //omega is the coefficient of the plaquette, divided by the module of the su2 submatrix for normalization
	  double omega_f=beta/(3*smod);
	  double z_norm=exp(-2*omega_f);
	  omega_f=1/omega_f;
	  
	  double temp_f,z_f,a0;
	  do
	    {
	      double z_temp=(z_norm-1)*rnd_get_unif(gen,0,1)+1;
	      a0     = 1+omega_f*log(z_temp);
	      z_f    = 1-a0*a0;
	      temp_f = sqr(rnd_get_unif(gen,0,1))-z_f;
	    }
	  while(temp_f>0);
	  
	  double x_rat=sqrt(z_f);
	  
	  //generate an su2 matrix
	  double fi=rnd_get_unif(gen,0,2*M_PI);
	  double cteta=rnd_get_unif(gen,-1,1);
	  double steta=sqrt(1-cteta*cteta);
	  
	  double a1=steta*cos(fi)*x_rat;
	  double a2=steta*sin(fi)*x_rat;
	  double a3=cteta*x_rat;
	  
	  double x0 = a0*r0 + a1*r1 + a2*r2 + a3*r3;
	  double x1 = r0*a1 - a0*r1 + a2*r3 - r2*a3;
	  double x2 = r0*a2 - a0*r2 + a3*r1 - r3*a1;
	  double x3 = r0*a3 - a0*r3 + a1*r2 - r1*a2;
	  
	  su2_prodassign_su3(x0,x1,x2,x3,isub_gr,prod);
	  su2_prodassign_su3(x0,x1,x2,x3,isub_gr,out);
	}
  }
  
  //return a single link after the overrelaxation procedure
  void su3_find_overrelaxed(su3 out,su3 in,su3 staple,int nov_hits)
  {
    //compute the original contribution to the action due to the given link 
    su3 prod;
    unsafe_su3_prod_su3_dag(prod,in,staple);
    
    //copy in link to out
    su3_copy(out,in);
    
    //iterate over overrelax hits
    for(int ihit=0;ihit<nov_hits;ihit++)
      //scan all the three possible subgroups
      for(int isub_gr=0;isub_gr<3;isub_gr++)
	{
	  //take the part of the su3 matrix
	  double r0,r1,r2,r3;
	  su2_part_of_su3(r0,r1,r2,r3,prod,isub_gr);
	  
	  //build the changing matrix
	  double x0=2*r0*r0-1;
	  double x1=-2*r0*r1;
	  double x2=-2*r0*r2;
	  double x3=-2*r0*r3;
	  
	  //change the link and optate the product
	  su2_prodassign_su3(x0,x1,x2,x3,isub_gr,prod);
	  su2_prodassign_su3(x0,x1,x2,x3,isub_gr,out);
	}
  }
  
  //return a cooled copy of the passed link
  void su3_find_cooled(su3 u,quad_su3 **eo_conf,int par,int ieo,int mu)
  {
    //compute the staple
    su3 staple;
    compute_point_summed_squared_staples_eo_conf_single_dir(staple,eo_conf,loclx_of_loceo[par][ieo],mu);
    
    //find the link that maximize the plaquette
    su3_unitarize_orthonormalizing(u,staple);
  }
  
  ////////////////////// products between su3 and color //////////////////
  
  //product of an su3 matrix by a color vector
  void unsafe_su3_prod_color(color a,su3 b,color c)
  {
    for(int c1=0;c1<3;c1++)
      {
	unsafe_complex_prod(a[c1],b[c1][0],c[0]);
	for(int c2=1;c2<3;c2++) complex_summ_the_prod(a[c1],b[c1][c2],c[c2]);
      }
  }
  
  //safe prod
  void safe_su3_prod_color(color a,su3 b,color c) {color t;unsafe_su3_prod_color(t,b,c);color_copy(a,t);}
  
  //summ
  void su3_summ_the_prod_color(color a,su3 b,color c)
  {
    for(int c1=0;c1<3;c1++)
      for(int c2=0;c2<3;c2++)
	complex_summ_the_prod(a[c1],b[c1][c2],c[c2]);
  }
  
  //subt
  void su3_subt_the_prod_color(color a,su3 b,color c)
  {
    for(int c1=0;c1<3;c1++)
      for(int c2=0;c2<3;c2++)
	complex_subt_the_prod(a[c1],b[c1][c2],c[c2]);
  }
  
  //dag prod
  void unsafe_su3_dag_prod_color(color a,su3 b,color c)
  {
    for(int c1=0;c1<3;c1++)
      {
	unsafe_complex_conj1_prod(a[c1],b[0][c1],c[0]);
	for(int c2=1;c2<3;c2++) complex_summ_the_conj1_prod(a[c1],b[c2][c1],c[c2]);
      }
  }
  
  //safe dag prod
  void safe_su3_dag_prod_color(color a,su3 b,color c)
  {color t;unsafe_su3_dag_prod_color(t,b,c);color_copy(a,t);}
  
  //summ dag
  void su3_dag_summ_the_prod_color(color a,su3 b,color c)
  {
    for(int c1=0;c1<3;c1++)
      for(int c2=0;c2<3;c2++) 
	complex_summ_the_conj1_prod(a[c1],b[c2][c1],c[c2]);
  }
  
  //subt dag
  void su3_dag_subt_the_prod_color(color a,su3 b,color c)
  {
    for(int c1=0;c1<3;c1++)
      for(int c2=0;c2<3;c2++)
	complex_subt_the_conj1_prod(a[c1],b[c2][c1],c[c2]);
  }
  
  //////////////////////////////////////////// color prod su3 ///////////////////////////////////////////
  
  //prod
  void unsafe_color_prod_su3(color a,color b,su3 c)
  {
    for(int c1=0;c1<3;c1++)
      {
	unsafe_complex_prod(a[c1],b[0],c[0][c1]);
	for(int c2=1;c2<3;c2++) complex_summ_the_prod(a[c1],b[c2],c[c2][c1]);
      }
  }
  
  //safe prod
  void safe_color_prod_su3(color a,color b,su3 c)
  {color t;unsafe_color_prod_su3(t,b,c);color_copy(a,t);}
  
  //dag
  void unsafe_color_prod_su3_dag(color a,color b,su3 c)
  {
    for(int c1=0;c1<3;c1++)
      {
	unsafe_complex_conj2_prod(a[c1],b[0],c[c1][0]);
	for(int c2=1;c2<3;c2++) complex_summ_the_conj2_prod(a[c1],b[c2],c[c1][c2]);
      }
  }
  
  //////////////////////////////////// Color and complex //////////////////////////////
  
  void safe_color_prod_complex(color out,color in,complex factor)
  {for(int i=0;i<3;i++) safe_complex_prod(((complex*)out)[i],((complex*)in)[i],factor);}
  void unsafe_color_prod_complex(color out,color in,complex factor)
  {for(int i=0;i<3;i++) unsafe_complex_prod(((complex*)out)[i],((complex*)in)[i],factor);}
  
  void safe_color_prod_complex_conj(color out,color in,complex factor)
  {for(int i=0;i<3;i++) safe_complex_conj2_prod(((complex*)out)[i],((complex*)in)[i],factor);}
  void unsafe_color_prod_complex_conj(color out,color in,complex factor)
  {for(int i=0;i<3;i++) unsafe_complex_conj2_prod(((complex*)out)[i],((complex*)in)[i],factor);}
  
  ////////////////////////////////// Operations between spincolor ////////////////////////
  
  //just print a spincolor
  void spincolor_print(spincolor c)
  {
    for(int id=0;id<4;id++)
      {
	for(int ic=0;ic<3;ic++) printf("%+016.16le,%+016.16le\t",c[id][ic][0],c[id][ic][1]);
	master_printf("\n");
      }
    master_printf("\n");
  }
  
  //summ two spincolors
  void spincolor_summ(spincolor a,spincolor b,spincolor c)
  {for(int i=0;i<24;i++) ((double*)a)[i]=((double*)b)[i]+((double*)c)[i];}
  void spincolor_summassign(spincolor a,spincolor b) {spincolor_summ(a,a,b);}
  
  //subtract two spincolors
  void spincolor_subt(spincolor a,spincolor b,spincolor c)
  {for(int i=0;i<24;i++) ((double*)a)[i]=((double*)b)[i]-((double*)c)[i];}
  void spincolor_subtassign(spincolor a,spincolor b) {spincolor_subt(a,a,b);}
  
  //summ after multyplying for a complex factor
  void safe_spincolor_summ_with_cfactor(spincolor a,spincolor b,spincolor c,complex factor)
  {
    complex temp;
    for(int i=0;i<12;i++)
      {
	unsafe_complex_prod(temp,factor,((complex*)c)[i]);
	complex_summ(((complex*)a)[i],((complex*)b)[i],temp);
      }
  }
  
  //summ after multyplying for a real factor
  void spincolor_summ_the_prod_double(spincolor a,spincolor b,spincolor c,double factor)
  {for(int i=0;i<24;i++) ((double*)a)[i]=factor*((double*)b)[i]+((double*)c)[i];}
  
  //spincolor*real
  void spincolor_prod_double(spincolor out,spincolor in,double factor)
  {for(int i=0;i<24;i++) ((double*)out)[i]=((double*)in)[i]*factor;}
  
  //spincolor*complex
  void unsafe_spincolor_prod_complex(spincolor out,spincolor in,complex factor)
  {for(int i=0;i<12;i++) unsafe_complex_prod(((complex*)out)[i],((complex*)in)[i],factor);}
  void safe_spincolor_prod_complex(spincolor out,spincolor in,complex factor)
  {for(int i=0;i<12;i++) safe_complex_prod(((complex*)out)[i],((complex*)in)[i],factor);}
  
  //spincolor+spincolor*complex
  void spincolor_summ_the_prod_complex(spincolor out,spincolor in,complex factor)
  {for(int i=0;i<12;i++) complex_summ_the_prod(((complex*)out)[i],((complex*)in)[i],factor);}
  
  //spincolor*i*real
  void unsafe_spincolor_summassign_the_prod_idouble(spincolor out,spincolor in,double factor)
  {
    for(int i=0;i<24;i+=2)
      {
	((double*)out)[i  ]-=((double*)in)[i+1]*factor;
	((double*)out)[i+1]+=((double*)in)[i  ]*factor;
      }
  }
  void spincolor_prodassign_idouble(spincolor out,double factor)
  {
    for(int i=0;i<24;i+=2)
      {
	double temp=-((double*)out)[i+1]*factor;
	((double*)out)[i+1]=((double*)out)[i  ]*factor;
	((double*)out)[i  ]=temp;
      }
  }
  
  //spincolor*i*real
  void unsafe_spincolor_summ_with_ifactor(spincolor out,spincolor a,spincolor b,double factor)
  {
    for(int i=0;i<24;i+=2)
      {
	((double*)out)[i  ]=((double*)a)[i  ]-((double*)b)[i+1]*factor;
	((double*)out)[i+1]=((double*)a)[i+1]+((double*)b)[i  ]*factor;
      }
  }
  
  //dirac*spincolor
  void unsafe_dirac_prod_spincolor(spincolor out,dirac_matr *m,spincolor in)
  {for(int id1=0;id1<4;id1++) for(int ic=0;ic<3;ic++) safe_complex_prod(out[id1][ic],m->entr[id1],in[m->pos[id1]][ic]);}
  
  //spincolor*dirac
  void unsafe_spincolor_prod_dirac(spincolor out,spincolor in,dirac_matr *m)
  {unsafe_dirac_prod_spincolor(out,m,in);}
  
  //dirac*spincolor
  void safe_dirac_prod_spincolor(spincolor out,dirac_matr *m,spincolor in)
  {
    spincolor tmp;
    unsafe_dirac_prod_spincolor(tmp,m,in);
    spincolor_copy(out,tmp);
  }
  
  //spincolor*dirac
  void safe_spincolor_prod_dirac(spincolor out,spincolor in,dirac_matr *m)
  {
    spincolor tmp;
    unsafe_spincolor_prod_dirac(tmp,in,m);
    spincolor_copy(out,tmp);
  }
  
  //su3*spincolor
  void unsafe_su3_prod_spincolor(spincolor out,su3 U,spincolor in)
  {for(int is=0;is<4;is++) unsafe_su3_prod_color(out[is],U,in[is]);}
  void su3_summ_the_prod_spincolor(spincolor out,su3 U,spincolor in)
  {for(int is=0;is<4;is++) su3_summ_the_prod_color(out[is],U,in[is]);}
  void su3_subt_the_prod_spincolor(spincolor out,su3 U,spincolor in)
  {for(int is=0;is<4;is++) su3_subt_the_prod_color(out[is],U,in[is]);}
  
  //su3^*spincolor
  void unsafe_su3_dag_prod_spincolor(spincolor out,su3 U,spincolor in)
  {for(int is=0;is<4;is++) unsafe_su3_dag_prod_color(out[is],U,in[is]);}
  void safe_su3_dag_prod_spincolor(spincolor out,su3 U,spincolor in)
  {spincolor temp;unsafe_su3_dag_prod_spincolor(temp,U,in);spincolor_copy(out,temp);}
  
  void su3_dag_summ_the_prod_spincolor(spincolor out,su3 U,spincolor in)
  {for(int is=0;is<4;is++) su3_dag_summ_the_prod_color(out[is],U,in[is]);}
  void su3_dag_subt_the_prod_spincolor(spincolor out,su3 U,spincolor in)
  {for(int is=0;is<4;is++) su3_dag_subt_the_prod_color(out[is],U,in[is]);}
  
  //su3^*gamma*spincolor
  void unsafe_su3_dag_dirac_prod_spincolor(spincolor out,su3 U,dirac_matr *m,spincolor in)
  {
    color tmp;
    for(int id1=0;id1<4;id1++)
      {
	for(int ic=0;ic<3;ic++) unsafe_complex_prod(tmp[ic],m->entr[id1],in[m->pos[id1]][ic]);
	unsafe_su3_dag_prod_color(out[id1],U,tmp);
      }
  }
  
  void unsafe_su3_dag_dirac_summ_the_prod_spincolor(spincolor out,su3 U,dirac_matr *m,spincolor in)
  {
    color tmp;
    for(int id1=0;id1<4;id1++)
      {
	for(int ic=0;ic<3;ic++) unsafe_complex_prod(tmp[ic],m->entr[id1],in[m->pos[id1]][ic]);
	su3_dag_summ_the_prod_color(out[id1],U,tmp);
      }
  }
  
  //su3*dirac*spincolor
  void unsafe_su3_dirac_prod_spincolor(spincolor out,su3 U,dirac_matr *m,spincolor in)
  {
    color tmp;
    for(int id1=0;id1<4;id1++)
      {
	for(int ic=0;ic<3;ic++) unsafe_complex_prod(tmp[ic],m->entr[id1],in[m->pos[id1]][ic]);
	unsafe_su3_prod_color(out[id1],U,tmp);
      }
  }
  
  void unsafe_su3_dirac_subt_the_prod_spincolor(spincolor out,su3 U,dirac_matr *m,spincolor in)
  {
    color tmp;
    for(int id1=0;id1<4;id1++)
      {
	for(int ic=0;ic<3;ic++) unsafe_complex_prod(tmp[ic],m->entr[id1],in[m->pos[id1]][ic]);
	su3_subt_the_prod_color(out[id1],U,tmp);
      }
  }
  
  ///////////////////////////////// su3*colorspinspin ///////////////////////////////////
  
  //colorspinspin*real
  void colorspinspin_prod_double(colorspinspin out,colorspinspin in,double factor)
  {for(int i=0;i<96;i++) ((double*)out)[i]=((double*)in)[i]*factor;}
  
  void unsafe_su3_prod_colorspinspin(colorspinspin a,su3 b,colorspinspin c)
  {
    for(int id_so=0;id_so<4;id_so++)
      for(int id_si=0;id_si<4;id_si++)
	for(int c1=0;c1<3;c1++)
	  {
	    unsafe_complex_prod(a[c1][id_si][id_so],b[c1][0],c[0][id_si][id_so]);
	    for(int c2=1;c2<3;c2++) complex_summ_the_prod(a[c1][id_si][id_so],b[c1][c2],c[c2][id_si][id_so]);
	  }
  }
  
  void unsafe_su3_dag_prod_colorspinspin(colorspinspin a,su3 b,colorspinspin c)
  {
    for(int id_so=0;id_so<4;id_so++)
      for(int id_si=0;id_si<4;id_si++)
	for(int c1=0;c1<3;c1++)
	  {
	    unsafe_complex_conj1_prod(a[c1][id_si][id_so],b[0][c1],c[0][id_si][id_so]);
	    for(int c2=1;c2<3;c2++)
	      complex_summ_the_conj1_prod(a[c1][id_si][id_so],b[c2][c1],c[c2][id_si][id_so]);
	  }
  }
  
  void su3_summ_the_prod_colorspinspin(colorspinspin a,su3 b,colorspinspin c)
  {
    for(int id_so=0;id_so<4;id_so++)
      for(int id_si=0;id_si<4;id_si++)
	for(int c1=0;c1<3;c1++)
	  for(int c2=0;c2<3;c2++)
	    complex_summ_the_prod(a[c1][id_si][id_so],b[c1][c2],c[c2][id_si][id_so]);
  }
  
  void su3_dag_subt_the_prod_colorspinspin(colorspinspin a,su3 b,colorspinspin c)
  {
    for(int id_so=0;id_so<4;id_so++)
      for(int id_si=0;id_si<4;id_si++)
	for(int c1=0;c1<3;c1++)
	  for(int c2=0;c2<3;c2++)
	    complex_subt_the_conj1_prod(a[c1][id_si][id_so],b[c2][c1],c[c2][id_si][id_so]);
  }
  
  void su3_dag_summ_the_prod_colorspinspin(colorspinspin a,su3 b,colorspinspin c)
  {
    for(int id_so=0;id_so<4;id_so++)
      for(int id_si=0;id_si<4;id_si++)
	for(int c1=0;c1<3;c1++)
	  for(int c2=0;c2<3;c2++)
	    complex_summ_the_conj1_prod(a[c1][id_si][id_so],b[c2][c1],c[c2][id_si][id_so]);
  }
  
  ////////////////////////////////////// su3spinspin /////////////////////////////////////
  
  //su3spinspin*complex
  void unsafe_su3spinspin_prod_complex(su3spinspin out,su3spinspin in,complex factor)
  {for(int i=0;i<144;i++) unsafe_complex_prod(((complex*)out)[i],((complex*)in)[i],factor);}
  void safe_su3spinspin_prod_complex(su3spinspin out,su3spinspin in,complex factor)
  {for(int i=0;i<144;i++) safe_complex_prod(((complex*)out)[i],((complex*)in)[i],factor);}
  
  void su3spinspin_prod_double(su3spinspin out,su3spinspin in,double factor)
  {for(int i=0;i<288;i++) ((double*)out)[i]=((double*)in)[i]*factor;}
  
  void unsafe_su3_prod_su3spinspin(su3spinspin a,su3 b,su3spinspin c)
  {
    for(int id_so=0;id_so<4;id_so++)
      for(int id_si=0;id_si<4;id_si++)
	for(int c1=0;c1<3;c1++)
	  for(int c3=0;c3<3;c3++)
	    {
	      unsafe_complex_prod(a[c1][c3][id_si][id_so],b[c1][0],c[0][c3][id_si][id_so]);
	      for(int c2=1;c2<3;c2++) complex_summ_the_prod(a[c1][c3][id_si][id_so],b[c1][c2],c[c2][c3][id_si][id_so]);
	    }
  }
  
  void unsafe_su3_dag_prod_su3spinspin(su3spinspin a,su3 b,su3spinspin c)
  {
    for(int id_so=0;id_so<4;id_so++)
      for(int id_si=0;id_si<4;id_si++)
	for(int c1=0;c1<3;c1++)
	  for(int c3=0;c3<3;c3++)
	    {
	      unsafe_complex_conj1_prod(a[c1][c3][id_si][id_so],b[0][c1],c[0][c3][id_si][id_so]);
	      for(int c2=1;c2<3;c2++)
		complex_summ_the_conj1_prod(a[c1][c3][id_si][id_so],b[c2][c1],c[c2][c3][id_si][id_so]);
	    }
  }
  
  void su3_summ_the_prod_su3spinspin(su3spinspin a,su3 b,su3spinspin c)
  {
    for(int id_so=0;id_so<4;id_so++)
      for(int id_si=0;id_si<4;id_si++)
	for(int c1=0;c1<3;c1++)
	  for(int c3=0;c3<3;c3++)
	    for(int c2=0;c2<3;c2++)
	      complex_summ_the_prod(a[c1][c3][id_si][id_so],b[c1][c2],c[c2][c3][id_si][id_so]);
  }
  
  void su3_dag_subt_the_prod_su3spinspin(su3spinspin a,su3 b,su3spinspin c)
  {
    for(int id_so=0;id_so<4;id_so++)
      for(int id_si=0;id_si<4;id_si++)
	for(int c1=0;c1<3;c1++)
	  for(int c2=0;c2<3;c2++)
	    for(int c3=0;c3<3;c3++)
	      complex_subt_the_conj1_prod(a[c1][c3][id_si][id_so],b[c2][c1],c[c2][c3][id_si][id_so]);
  }
  
  void su3_dag_summ_the_prod_su3spinspin(su3spinspin a,su3 b,su3spinspin c)
  {
    for(int id_so=0;id_so<4;id_so++)
      for(int id_si=0;id_si<4;id_si++)
	for(int c1=0;c1<3;c1++)
	  for(int c2=0;c2<3;c2++)
	    for(int c3=0;c3<3;c3++)
	      complex_summ_the_conj1_prod(a[c1][c3][id_si][id_so],b[c2][c1],c[c2][c3][id_si][id_so]);
  }
  
  //put a matrix to random used passed random generator
  void su3_put_to_rnd(su3 u_ran,rnd_gen &rnd)
  {
    su3_put_to_id(u_ran);
    
    for(int i1=0;i1<3;i1++)
      for(int i2=i1+1;i2<3;i2++)
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
