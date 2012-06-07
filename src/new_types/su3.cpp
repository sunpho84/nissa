#include <math.h>
#include <string.h>

#include "new_types_definitions.h"
#include "../base/routines.h"
#include "../base/random.h"
#include "float128.h"
#include "complex.h"
#include "su3.h"

//////////////////////////////////// Put to zero or 1 /////////////////////////////////

void color_put_to_zero(color m)
{memset(m,0,sizeof(color));}
void su3_put_to_zero(su3 m)
{memset(m,0,sizeof(su3));}
void as2t_su3_put_to_zero(as2t_su3 m)
{memset(m,0,sizeof(as2t_su3));}
void spincolor_put_to_zero(spincolor m)
{memset(m,0,sizeof(spincolor));}
void su3spinspin_put_to_zero(su3spinspin m)
{memset(m,0,sizeof(su3spinspin));}
void su3_put_to_id(su3 m)
{su3_put_to_zero(m);for(int ic=0;ic<3;ic++) m[ic][ic][0]=1;}

//////////////////////////////////////// Copy /////////////////////////////////////

void color_copy(color b,color a)
{memcpy(b,a,sizeof(color));}
void su3_copy(su3 b,su3 a)
{memcpy(b,a,sizeof(su3));}
void quad_su3_copy(quad_su3 b,quad_su3 a)
{memcpy(b,a,sizeof(quad_su3));}
void spincolor_copy(spincolor b,spincolor a)
{memcpy(b,a,sizeof(spincolor));}

////////////////////////////////// Operations between colors //////////////////////////

//just print a color
void color_print(color c)
{
  for(int icol=0;icol<3;icol++)
    {
      for(int ri=0;ri<2;ri++)
	{
	  master_printf("%16.16lg",c[icol][ri]);
	  if(ri==0) printf("\t");
	}
      master_printf("\n");
    }
  master_printf("\n");
}

//summ two colors
void color_summ(color a,color b,color c)
{for(int i=0;i<6;i++) ((double*)a)[i]=((double*)b)[i]+((double*)c)[i];}

void color_isumm(color a,color b,color c)
{for(int i=0;i<6;i+=2) {((double*)a)[i]=((double*)b)[i]-((double*)c)[i+1];((double*)a)[i+1]=((double*)b)[i+1]+((double*)c)[i];}}

void color_isubt(color a,color b,color c)
{for(int i=0;i<6;i+=2) {((double*)a)[i]=((double*)b)[i]+((double*)c)[i+1];((double*)a)[i+1]=((double*)b)[i+1]-((double*)c)[i];}}

void color_subt(color a,color b,color c)
{for(int i=0;i<6;i++) ((double*)a)[i]=((double*)b)[i]-((double*)c)[i];}

void color_summassign(color a,color b)
{for(int i=0;i<6;i++) ((double*)a)[i]+=((double*)b)[i];}

void color_subtassign(color a,color b)
{for(int i=0;i<6;i++) ((double*)a)[i]-=((double*)b)[i];}

void color_isummassign(color a,color b)
{for(int i=0;i<6;i+=2) {((double*)a)[i]-=((double*)b)[i+1];((double*)a)[i+1]+=((double*)b)[i];}}

void color_isubtassign(color a,color b)
{for(int i=0;i<6;i+=2) {((double*)a)[i]+=((double*)b)[i+1];((double*)a)[i+1]-=((double*)b)[i];}}

void color_prod_double(color a,color b,double c)
{for(int i=0;i<6;i++) ((double*)a)[i]=((double*)b)[i]*c;}

/////////////////////////////// Generate an hermitean matrix ///////////////////////

//Taken from M.D'elia
void herm_put_to_gauss(su3 H,rnd_gen *gen,double sigma)
{
  const double one_by_sqrt2=0.707106781186547;
  const double one_by_sqrt3=0.577350269189626;
  const double two_by_sqrt3=1.15470053837925;
  
  double r[8];
  for(int ir=0;ir<8;ir++)
    r[ir]=rnd_get_gauss(gen,0,sigma*one_by_sqrt2);

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

//just print an su3 matrix
void su3_print(su3 U)
{
  for(int icol1=0;icol1<3;icol1++)
    {
      for(int icol2=0;icol2<3;icol2++)
	{
	  for(int ri=0;ri<2;ri++)
	    {
	      master_printf("%16.16lg",U[icol1][icol2][ri]);
	      if(ri==0) printf("\t");
	    }
	  master_printf("\t");
	}
      master_printf("\n");
    }
  master_printf("\n");
}

//return the trace of an su3 matrix
void su3_trace(complex tr,su3 m)
{
  tr[0]=m[0][0][0];
  tr[1]=m[0][0][1];

  complex_summ(tr,tr,m[1][1]);
  complex_summ(tr,tr,m[2][2]);
}

//return the anti-hermitian traceless part of an su3 matrix
void su3_traceless_anti_hermitian_part(su3 out,su3 in)
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

//return the hemitian su3 matrix
void safe_su3_hermitian(su3 out,su3 in)
{
  su3 tmp;
  for(int ic_in=0;ic_in<3;ic_in++)
    for(int ic_out=0;ic_out<3;ic_out++)
      {
	tmp[ic_in][ic_out][0]= in[ic_out][ic_in][0];
	tmp[ic_in][ic_out][1]=-in[ic_out][ic_in][1];
      }
  su3_copy(out,tmp);
}
void unsafe_su3_hermitian(su3 out,su3 in)
{
  for(int ic_in=0;ic_in<3;ic_in++)
    for(int ic_out=0;ic_out<3;ic_out++)
      {
	out[ic_in][ic_out][0]= in[ic_out][ic_in][0];
	out[ic_in][ic_out][1]=-in[ic_out][ic_in][1];
      }
}

//summ two su3 matrixes
void su3_summ(su3 a,su3 b,su3 c)
{for(int i=0;i<18;i++) ((double*)a)[i]=((double*)b)[i]+((double*)c)[i];}
void su3_summ_real(su3 a,su3 b,double c)
{su3_copy(a,b);for(int i=0;i<3;i++) a[i][i][0]=b[i][i][0]+c;}
void su3_subt(su3 a,su3 b,su3 c)
{for(int i=0;i<18;i++) ((double*)a)[i]=((double*)b)[i]-((double*)c)[i];}
void su3_subt_complex(su3 a,su3 b,complex c)
{for(int i=0;i<3;i++) for(int ri=0;ri<2;ri++) a[i][i][ri]=b[i][i][ri]-c[ri];}

//Product of two su3 matrixes
void unsafe_su3_prod_su3(su3 a,su3 b,su3 c)
{
  for(int ir_out=0;ir_out<3;ir_out++)
    for(int ic_out=0;ic_out<3;ic_out++)
      {
	unsafe_complex_prod(a[ir_out][ic_out],b[ir_out][0],c[0][ic_out]);
	for(int itemp=1;itemp<3;itemp++)
	  complex_summ_the_prod(a[ir_out][ic_out],b[ir_out][itemp],c[itemp][ic_out]);
      }
}
void safe_su3_prod_su3(su3 a,su3 b,su3 c)
{su3 d;unsafe_su3_prod_su3(d,b,c);su3_copy(a,d);}

//Product of two su3 matrixes
void unsafe_su3_dag_prod_su3(su3 a,su3 b,su3 c)
{
  for(int ir_out=0;ir_out<3;ir_out++)
    for(int ic_out=0;ic_out<3;ic_out++)
      {
	unsafe_complex_conj1_prod(a[ir_out][ic_out],b[0][ir_out],c[0][ic_out]);
	for(int itemp=1;itemp<3;itemp++)
	  complex_summ_the_conj1_prod(a[ir_out][ic_out],b[itemp][ir_out],c[itemp][ic_out]);
      }
}
void safe_su3_dag_prod_su3(su3 a,su3 b,su3 c)
{su3 d;unsafe_su3_dag_prod_su3(d,b,c);su3_copy(a,d);}

//Product of two su3 matrixes
void unsafe_su3_prod_su3_dag(su3 a,su3 b,su3 c)
{
  for(int ir_out=0;ir_out<3;ir_out++)
    for(int ic_out=0;ic_out<3;ic_out++)
      {
	unsafe_complex_conj2_prod(a[ir_out][ic_out],b[ir_out][0],c[ic_out][0]);
	complex_summ_the_conj2_prod(a[ir_out][ic_out],b[ir_out][1],c[ic_out][1]);
	complex_summ_the_conj2_prod(a[ir_out][ic_out],b[ir_out][2],c[ic_out][2]);
      }
}
void safe_su3_prod_su3_dag(su3 a,su3 b,su3 c)
{su3 d;unsafe_su3_prod_su3_dag(d,b,c);su3_copy(a,d);}

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

//Product of two su3 matrixes
void unsafe_su3_dag_prod_su3_dag(su3 a,su3 b,su3 c)
{
  for(int ir_out=0;ir_out<3;ir_out++)
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
void su3_hermitian_prod_double(su3 a,su3 b,double r)
{
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      {
	a[i][j][0]= r*b[j][i][0];
	a[i][j][1]=-r*b[j][i][1];
      }
}

//summ the prod of su3 with imag
void su3_prod_with_idouble(su3 a,su3 b,double r)
{
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      {
	double c=b[i][j][0]*r;
	a[i][j][0]=-b[i][j][1]*r;
	a[i][j][1]=c;
      }
}

//summ the prod of su3 with real
void su3_summ_the_prod_double(su3 a,su3 b,double r)
{
  double *da=(double*)a;
  double *db=(double*)b;

  for(int i=0;i<18;i++) da[i]+=db[i]*r;
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
  //two indices specifying subgroup
  int sub_gr_index[3][2]={{0,1},{1,2},{0,2}};
  //loop over the three subgroups
  for(int isub_gr=0;isub_gr<3;isub_gr++)
    {
      int a=sub_gr_index[isub_gr][0];
      int b=sub_gr_index[isub_gr][1];
      
      //compute the product
      su3 P;
      unsafe_su3_prod_su3_dag(P,U,M);
      
      //take projection of prod
      double A=P[a][a][0]+P[b][b][0];
      double B=P[a][b][1]+P[b][a][1];
      double C=P[a][b][0]-P[b][a][0];
      double D=P[a][a][1]-P[b][b][1];
      
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
      
      //define the su2 matrix
      complex r[2][2]={{{A,-D},{-C,-B}},{{C,-B},{ A, D}}};
      
      //update the guess: one of the line do not change, the other two changes
      for(int c=0;c<3;c++)
	{
	  complex t1,t2;
	  
	  unsafe_complex_prod(  t1,r[0][0],U[a][c]);
	  complex_summ_the_prod(t1,r[0][1],U[b][c]);
	  
	  unsafe_complex_prod(  t2,r[1][0],U[a][c]);
	  complex_summ_the_prod(t2,r[1][1],U[b][c]);
	  
	  complex_copy(U[a][c],t1);
	  complex_copy(U[b][c],t2);
	}
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
void safe_su3_prod_color(color a,su3 b,color c)
{color t;unsafe_su3_prod_color(t,b,c);color_copy(a,t);}

//summ
void su3_summ_the_prod_color(color a,su3 b,color c)
{
  for(int c1=0;c1<3;c1++)
    {
      complex_summ_the_prod(a[c1],b[c1][0],c[0]);
      for(int c2=1;c2<3;c2++) complex_summ_the_prod(a[c1],b[c1][c2],c[c2]);
    }
}

//subt
void su3_subt_the_prod_color(color a,su3 b,color c)
{
  for(int c1=0;c1<3;c1++)
    {
      complex_subt_the_prod(a[c1],b[c1][0],c[0]);
      for(int c2=1;c2<3;c2++) complex_subt_the_prod(a[c1],b[c1][c2],c[c2]);
    }
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
    {
      complex_summ_the_conj1_prod(a[c1],b[0][c1],c[0]);
      for(int c2=1;c2<3;c2++) complex_summ_the_conj1_prod(a[c1],b[c2][c1],c[c2]);
    }
}

//subt dag
void su3_dag_subt_the_prod_color(color a,su3 b,color c)
{
  for(int c1=0;c1<3;c1++)
    {
      complex_subt_the_conj1_prod(a[c1],b[0][c1],c[0]);
      for(int c2=1;c2<3;c2++) complex_subt_the_conj1_prod(a[c1],b[c2][c1],c[c2]);
    }
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

////////////////////////////////// Operations between spincolor ///////////////////

//summ two spincolors
void spincolor_summ(spincolor a,spincolor b,spincolor c)
{for(int i=0;i<24;i++) ((double*)a)[i]=((double*)b)[i]+((double*)c)[i];}
void spincolor_summassign(spincolor a,spincolor b)
{spincolor_summ(a,a,b);}

//subtract two spincolors
void spincolor_subt(spincolor a,spincolor b,spincolor c)
{for(int i=0;i<24;i++) ((double*)a)[i]=((double*)b)[i]-((double*)c)[i];}
void spincolor_subtassign(spincolor a,spincolor b)
{spincolor_subt(a,a,b);}

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

void unsafe_su3_dag_summ_the_prod_spincolor(spincolor out,su3 U,spincolor in)
{for(int is=0;is<4;is++) su3_dag_summ_the_prod_color(out[is],U,in[is]);}
void unsafe_su3_dag_subt_the_prod_spincolor(spincolor out,su3 U,spincolor in)
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

//su3spinspin*complex
void unsafe_su3spinspin_prod_complex(su3spinspin out,su3spinspin in,complex factor)
{for(int i=0;i<144;i++) unsafe_complex_prod(((complex*)out)[i],((complex*)in)[i],factor);}
void safe_su3spinspin_prod_complex(su3spinspin out,su3spinspin in,complex factor)
{for(int i=0;i<144;i++) safe_complex_prod(((complex*)out)[i],((complex*)in)[i],factor);}
