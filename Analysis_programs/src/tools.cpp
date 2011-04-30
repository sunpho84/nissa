#pragma once

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

string combine(const char *format,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);
  
  return string(buffer);
}

//Estrae a real number between 0 and max
inline double rnd(double max=1)
{return random()*max/2147483648.0;}

//Extract a gaussian number
double rng(double med=0,double sigma=1)
{
  double q,r,x;
  static double y;
  static bool flag=true;

  if(flag)
    {
      r=sqrt(-2*log(1-rnd())); //This is correct, otherwise may be inf
      q=2*M_PI*rnd();
      
      x=r*cos(q);
      y=r*sin(q);
    }
  else x=y;

  flag=!flag;

  return med+sigma*x;
}
