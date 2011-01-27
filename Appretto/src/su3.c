#pragma once

#include "global.c"

//////////////////////////////////// Put to zero /////////////////////////////////

void color_put_to_zero(color m){memset(m,0,sizeof(color));}
void su3_put_to_zero(su3 m){memset(m,0,sizeof(su3));}
void as2t_su3_put_to_zero(as2t_su3 m){memset(m,0,sizeof(as2t_su3));}
void spincolor_put_to_zero(spincolor m){memset(m,0,sizeof(spincolor));}
void su3spinspin_put_to_zero(su3spinspin m){memset(m,0,sizeof(su3spinspin));}
void su3_put_to_id(su3 m){su3_put_to_zero(m);for(int ic=0;ic<3;ic++) m[ic][ic][0]=1;}

//////////////////////////////////////// Copy /////////////////////////////////////

void color_copy(color b,color a){memcpy(b,a,sizeof(color));}
void su3_copy(su3 b,su3 a){memcpy(b,a,sizeof(su3));}
void quad_su3_copy(quad_su3 b,quad_su3 a){memcpy(b,a,sizeof(quad_su3));}
void spincolor_copy(spincolor b,spincolor a){memcpy(b,a,sizeof(spincolor));}

////////////////////////////////// Operations between colors //////////////////////////

//summ two colors
void color_summ(color a,color b,color c)
{for(int i=0;i<6;i++) ((double*)a)[i]=((double*)b)[i]+((double*)c)[i];}

void color_isumm(color a,color b,color c)
{for(int i=0;i<6;i+=2) {((double*)a)[i]=((double*)b)[i]-((double*)c)[i+1];((double*)a)[i+1]=((double*)b)[i+1]+((double*)c)[i];}}

void color_isubt(color a,color b,color c)
{for(int i=0;i<6;i+=2) {((double*)a)[i]=((double*)b)[i]+((double*)c)[i+1];((double*)a)[i+1]=((double*)b)[i+1]-((double*)c)[i];}}

void color_subt(color a,color b,color c)
{for(int i=0;i<6;i++) ((double*)a)[i]=((double*)b)[i]-((double*)c)[i];}

void summassign_color(color a,color b)
{for(int i=0;i<6;i++) ((double*)a)[i]+=((double*)b)[i];}

void subtassign_color(color a,color b)
{for(int i=0;i<6;i++) ((double*)a)[i]-=((double*)b)[i];}

void summassign_icolor(color a,color b)
{for(int i=0;i<6;i+=2) {((double*)a)[i]-=((double*)b)[i+1];((double*)a)[i+1]+=((double*)b)[i];}}

void subtassign_icolor(color a,color b)
{for(int i=0;i<6;i+=2) {((double*)a)[i]+=((double*)b)[i+1];((double*)a)[i+1]-=((double*)b)[i];}}

////////////////////////////////// Operations between su3 //////////////////////////

//return the trace of an su3 matrix
void su3_trace(complex tr,su3 m)
{
  tr[0]=m[0][0][0];
  tr[1]=m[0][0][1];

  complex_summ(tr,tr,m[1][1]);
  complex_summ(tr,tr,m[2][2]);
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

//Product of two su3 matrixes
void su3_prod_su3(su3 a,su3 b,su3 c)
{
  for(int ir_out=0;ir_out<3;ir_out++)
    for(int ic_out=0;ic_out<3;ic_out++)
      {
	unsafe_complex_prod(a[ir_out][ic_out],b[ir_out][0],c[0][ic_out]);
	for(int itemp=1;itemp<3;itemp++)
	  complex_summ_the_prod(a[ir_out][ic_out],b[ir_out][itemp],c[itemp][ic_out]);
      }
}

//Product of two su3 matrixes
void su3_dag_prod_su3(su3 a,su3 b,su3 c)
{
  for(int ir_out=0;ir_out<3;ir_out++)
    for(int ic_out=0;ic_out<3;ic_out++)
      {
	unsafe_complex_conj1_prod(a[ir_out][ic_out],b[0][ir_out],c[0][ic_out]);
	for(int itemp=1;itemp<3;itemp++)
	  complex_summ_the_conj1_prod(a[ir_out][ic_out],b[itemp][ir_out],c[itemp][ic_out]);
      }
}

//Product of two su3 matrixes
void su3_prod_su3_dag(su3 a,su3 b,su3 c)
{
  for(int ir_out=0;ir_out<3;ir_out++)
    for(int ic_out=0;ic_out<3;ic_out++)
      {
	unsafe_complex_conj2_prod(a[ir_out][ic_out],b[ir_out][0],c[ic_out][0]);
	for(int itemp=1;itemp<3;itemp++)
	  complex_summ_the_conj2_prod(a[ir_out][ic_out],b[ir_out][itemp],c[ic_out][itemp]);
      }
}

//Product of two su3 matrixes
void su3_dag_prod_su3_dag(su3 a,su3 b,su3 c)
{
  for(int ir_out=0;ir_out<3;ir_out++)
    for(int ic_out=0;ic_out<3;ic_out++)
      {
	unsafe_complex_conj_conj_prod(a[ir_out][ic_out],b[0][ir_out],c[ic_out][0]);
	for(int itemp=1;itemp<3;itemp++)
	  complex_summ_the_conj_conj_prod(a[ir_out][ic_out],b[itemp][ir_out],c[ic_out][itemp]);
      }
}

//product of an su3 matrix by a complex
void safe_su3_prod_complex(su3 a,su3 b,complex c)
{
  complex *ca=(complex*)a;
  complex *cb=(complex*)b;

  for(int i=0;i<9;i++) safe_complex_prod(ca[i],cb[i],c);
}

//product of an su3 matrix by a real
void su3_prod_real(su3 a,su3 b,double r)
{
  double *da=(double*)a;
  double *db=(double*)b;

  for(int i=0;i<18;i++) da[i]=db[i]*r;
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

void unsafe_color_prod_su3(color a,color b,su3 c)
{
  for(int c1=0;c1<3;c1++)
    {
      unsafe_complex_prod(a[c1],b[0],c[0][c1]);
      for(int c2=1;c2<3;c2++) complex_summ_the_prod(a[c1],b[c2],c[c2][c1]);
    }
}

void unsafe_summ_su3_prod_color(color a,su3 b,color c)
{for(int c1=0;c1<3;c1++) for(int c2=0;c2<3;c2++) complex_summ_the_prod(a[c1],b[c1][c2],c[c2]);}

void unsafe_subt_su3_prod_color(color a,su3 b,color c)
{for(int c1=0;c1<3;c1++) for(int c2=0;c2<3;c2++) complex_subt_the_prod(a[c1],b[c1][c2],c[c2]);}

//product of an su3 matrix by a color vector
void safe_su3_prod_color(color a,su3 b,color c)
{
  color t;
  for(int c1=0;c1<3;c1++)
    {
      unsafe_complex_prod(t[c1],b[c1][0],c[0]);
      for(int c2=1;c2<3;c2++) complex_summ_the_prod(t[c1],b[c1][c2],c[c2]);
    }
  color_copy(a,t);
}

//product of an su3 matrix by a color vector
void unsafe_su3_dag_prod_color(color a,su3 b,color c)
{
  for(int c1=0;c1<3;c1++)
    {
      unsafe_complex_conj1_prod(a[c1],b[0][c1],c[0]);
      for(int c2=1;c2<3;c2++) complex_summ_the_conj1_prod(a[c1],b[c2][c1],c[c2]);
    }
}

void unsafe_color_prod_su3_dag(color a,color b,su3 c)
{
  for(int c1=0;c1<3;c1++)
    {
      unsafe_complex_conj2_prod(a[c1],b[0],c[c1][0]);
      for(int c2=1;c2<3;c2++) complex_summ_the_conj2_prod(a[c1],b[c2],c[c1][c2]);
    }
}

void unsafe_summ_su3_dag_prod_color(color a,su3 b,color c)
{for(int c1=0;c1<3;c1++) for(int c2=0;c2<3;c2++) complex_summ_the_conj1_prod(a[c1],b[c2][c1],c[c2]);}

//product of an su3 matrix by a color vector
void safe_su3_dag_prod_color(color a,su3 b,color c)
{
  color t;
  for(int c1=0;c1<3;c1++)
    {
      safe_complex_conj1_prod(t[c1],b[0][c1],c[0]);
      for(int c2=1;c2<3;c2++) complex_summ_the_conj1_prod(t[c1],b[c2][c1],c[c2]);
    }
  color_copy(a,t);
}

////////////////////////////////// Operations between spincolor ///////////////////

//summ two spincolors
void spincolor_summ(spincolor a,spincolor b,spincolor c)
{for(int i=0;i<24;i++) ((double*)a)[i]=((double*)b)[i]+((double*)c)[i];}

//subtract two spincolors
void spincolor_subt(spincolor a,spincolor b,spincolor c)
{for(int i=0;i<24;i++) ((double*)a)[i]=((double*)b)[i]-((double*)c)[i];}

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
void safe_spincolor_summ_with_rfactor(spincolor a,spincolor b,spincolor c,double factor)
{for(int i=0;i<24;i++) ((double*)a)[i]=factor*((double*)b)[i]+((double*)c)[i];}

//spincolor*complex
void safe_complex_prod_spincolor(spincolor out,complex a,spincolor in)
{for(int i=0;i<12;i++) safe_complex_prod(((complex*)out)[i],a,((complex*)in)[i]);}

//spincolor*real
void assign_spincolor_prod_real(spincolor out,double factor)
{for(int i=0;i<24;i++) ((double*)out)[i]*=factor;}

//summ assign
void summassign_spincolor(spincolor out,spincolor in)
{for(int i=0;i<24;i++) ((double*)out)[i]+=((double*)in)[i];}

//subt assign
void subtassign_spincolor(spincolor out,spincolor in)
{for(int i=0;i<24;i++) ((double*)out)[i]-=((double*)in)[i];}

//spincolor*complex
void unsafe_spincolor_prod_complex(spincolor out,spincolor in,complex factor)
{for(int i=0;i<12;i++) unsafe_complex_prod(((complex*)out)[i],((complex*)in)[i],factor);}

//spincolor*complex
void spincolor_summ_the_prod_complex(spincolor out,spincolor in,complex factor)
{for(int i=0;i<12;i++) complex_summ_the_prod(((complex*)out)[i],((complex*)in)[i],factor);}

//spincolor*real
void unsafe_summassign_spincolor_prod_real(spincolor out,spincolor in,double factor)
{for(int i=0;i<24;i++) ((double*)out)[i]+=((double*)in)[i]*factor;}

//spincolor*i*real
void unsafe_summassign_spincolor_prod_ireal(spincolor out,spincolor in,double factor)
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

//dirac*spincolor
void safe_dirac_prod_spincolor(spincolor out,dirac_matr *m,spincolor in)
{
  spincolor tmp;
  unsafe_dirac_prod_spincolor(tmp,m,in);
  spincolor_copy(out,tmp);
}

//su3*spincolor
void unsafe_su3_prod_spincolor(spincolor out,su3 U,spincolor in)
{for(int is=0;is<4;is++) unsafe_su3_prod_color(out[is],U,in[is]);}

void unsafe_summ_su3_prod_spincolor(spincolor out,su3 U,spincolor in)
{for(int is=0;is<4;is++) unsafe_summ_su3_prod_color(out[is],U,in[is]);}

//su3^*spincolor
void unsafe_su3_dag_prod_spincolor(spincolor out,su3 U,spincolor in)
{for(int is=0;is<4;is++) unsafe_su3_dag_prod_color(out[is],U,in[is]);}

void unsafe_summ_su3_dag_prod_spincolor(spincolor out,su3 U,spincolor in)
{for(int is=0;is<4;is++) unsafe_summ_su3_dag_prod_color(out[is],U,in[is]);}

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

void unsafe_summ_su3_dag_dirac_prod_spincolor(spincolor out,su3 U,dirac_matr *m,spincolor in)
{
  color tmp;
  for(int id1=0;id1<4;id1++)
    {
      for(int ic=0;ic<3;ic++) unsafe_complex_prod(tmp[ic],m->entr[id1],in[m->pos[id1]][ic]);
      unsafe_summ_su3_dag_prod_color(out[id1],U,tmp);
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

void unsafe_subt_su3_dirac_prod_spincolor(spincolor out,su3 U,dirac_matr *m,spincolor in)
{
  color tmp;
  for(int id1=0;id1<4;id1++)
    {
      for(int ic=0;ic<3;ic++) unsafe_complex_prod(tmp[ic],m->entr[id1],in[m->pos[id1]][ic]);
      unsafe_subt_su3_prod_color(out[id1],U,tmp);
    }
}
