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
  
  int i=0;
  
#ifdef BGP
  bgp_complex N0,N1,N2;
  bgp_color_put_to_zero(N0,N1,N2);
  while(i+24<n)
    {
      bgp_complex A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,AA,AB;
      bgp_complex B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,BA,BB;
      bgp_spincolor_load(A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,AA,AB,(*((spincolor*)(a+i))));
      bgp_spincolor_load(B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,BA,BB,(*((spincolor*)(b+i))));
      bgp_summassign_color_realscalarprod_spincolor(N0,N1,N2, A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,AA,AB, B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,BA,BB);
      
      i+=24;
    }
  bgp_square_norm_color(N0,N1,N2);
  complex cout;
  bgp_complex_save(cout,N0);
  out=cout[0];
#endif
  
  //reduce the remaining data (or all if not on bgp)
  while(i<n)
    {
      out+=a[i]*b[i];
      i++;
    }
  
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
