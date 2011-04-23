#pragma once

double unary_minus(const double a){return -a;}
double unary_plus(const double a){return a;}
double sqr(const double a){return a*a;}

double double_summ(const double a,const double b){return a+b;}
double double_subt(const double a,const double b){return a-b;}
double double_prod(const double a,const double b){return a*b;}
double double_frac(const double a,const double b){return a/b;}

/////////////////////////////////////////////////////////////

double lin_fun(double const x,double*p)
{return p[0]*x+p[1];}

double const_fun(double const x,double*p)
{return p[0];}
