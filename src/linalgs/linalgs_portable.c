#pragma once

//set to zero
void double_vector_init_to_zero(double *a,int n)
{
  memset(a,0,n*sizeof(double));
  set_borders_invalid(a);
}

//copy
void double_vector_copy(double *a,double *b,int n)
{
  memcpy(a,b,n*sizeof(double));
  set_borders_invalid(a);
}

//scalar product between a and b
double double_vector_loc_scalar_prod(double *a,double *b,int n)
{
  double out=0;
  
  for(int i=0;i<n;i++) out+=a[i]*b[i];
  
  return out;
}

//a[]=b[]+c[]*d
void double_vector_summ_double_vector_prod_double(double *a,double *b,double *c,double d,int n)
{
  for(int i=0;i<n;i++) a[i]=b[i]+c[i]*d;
  set_borders_invalid(a);
}

//a[]=b[]*c+d[]*e
void double_vector_linear_comb(double *a,double *b,double c,double *d,double e,int n)
{
  for(int i=0;i<n;i++) a[i]=b[i]*c+d[i]*e;
  set_borders_invalid(a);
}

//a[]=b[]-c[]*d
void double_vector_subt_double_vector_prod_double(double *a,double *b,double *c,double d,int n)
{double_vector_summ_double_vector_prod_double(a,b,c,-d,n);}

//global reduction wrapper
double double_vector_glb_scalar_prod(double *a,double *b,int n)
{return glb_reduce_double(double_vector_loc_scalar_prod(a,b,n));}
