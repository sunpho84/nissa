#pragma once

#define njack 10

typedef double jack[njack+1];

double jack_error(jack data)
{
  double sx=0,s2x=0;
  
  for(int ij=0;ij<njack;ij++)
    {
      sx+=data[ij];
      s2x+=data[ij]*data[ij];
    }
  sx/=njack;
  s2x/=njack;
  s2x-=sx*sx;
  
  return sqrt(s2x*(njack-1)); 
}

void jack_prod_jack(jack c,jack a,jack b)
{
  for(int ijack=0;ijack<njack+1;ijack++)
    c[ijack]=a[ijack]*b[ijack];
}

void jack_fracassign_jack(jack c,jack a)
{
  for(int ijack=0;ijack<njack+1;ijack++)
    c[ijack]/=a[ijack];
}

void jack_prodassign_jack(jack c,jack a)
{
  for(int ijack=0;ijack<njack+1;ijack++)
    c[ijack]*=a[ijack];
}

void jack_frac_jack(jack c,jack a,jack b)
{
  for(int ijack=0;ijack<njack+1;ijack++)
    c[ijack]=a[ijack]/b[ijack];
}

void double_frac_jack(jack c,double a,jack b)
{
  for(int ijack=0;ijack<njack+1;ijack++)
    c[ijack]=a/b[ijack];
}

void jack_subtprod_jack(jack c,jack a,jack b)
{
  for(int ijack=0;ijack<njack+1;ijack++)
    c[ijack]-=a[ijack]*b[ijack];
}

void jack_subt_jack(jack c,jack a,jack b)
{
  for(int ijack=0;ijack<njack+1;ijack++)
    c[ijack]=a[ijack]-b[ijack];
}

void jack_summ_jack(jack c,jack a,jack b)
{
  for(int ijack=0;ijack<njack+1;ijack++)
    c[ijack]=a[ijack]+b[ijack];
}

void jack_summassign_jack(jack c,jack a)
{
  for(int ijack=0;ijack<njack+1;ijack++)
    c[ijack]+=a[ijack];
}

void jack_fracassign_double(jack c,double a)
{
  for(int ijack=0;ijack<njack+1;ijack++)
    c[ijack]/=a;
}

void jack_prodassign_double(jack c,double a)
{
  for(int ijack=0;ijack<njack+1;ijack++)
    c[ijack]*=a;
}

void jack_summassign_double(jack c,double a)
{
  for(int ijack=0;ijack<njack+1;ijack++)
    c[ijack]+=a;
}

void jack_summ_double(jack c,jack a,double b)
{
  for(int ijack=0;ijack<njack+1;ijack++)
    c[ijack]=a[ijack]+b;
}

void jack_subt_double(jack c,jack a,double b)
{
  for(int ijack=0;ijack<njack+1;ijack++)
    c[ijack]=a[ijack]-b;
}

void jack_from_double(jack c,double b)
{
  for(int ijack=0;ijack<njack+1;ijack++)
    c[ijack]=b;
}

void jack_sqrt_jack(jack c,jack a)
{
  for(int ijack=0;ijack<njack+1;ijack++)
    c[ijack]=sqrt(a[ijack]);
}

