#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#if HIGH_PREC_TYPE==GMP_HIGH_PREC
 #include <gmpxx.h>
#endif

#include "base/debug.hpp"
#include "routines/ios.hpp"
#include "high_prec.hpp"

namespace nissa
{
  //initialize the high precision
  void init_high_precision()
  {
#if HIGH_PREC_TYPE==GMP_HIGH_PREC
    //init default precision for gmp
    mpf_set_default_prec(mpf_precision);
    master_printf("Support for >128 bit precision: GMP\n");
#else
    master_printf("Support for >128 bit precision: NATIVE\n");
#endif
    
    //perform a sanity check on float 128
    check_128_bit_prec();
  }
  
#if (HIGH_PREC_TYPE==NATIVE_HIGH_PREC)
  int high_prec_nbits() {return 209;}
  void float_high_prec_t_print(float_high_prec_t a){master_printf("%lg %lg %lg %lg\n",a[0],a[1],a[2],a[3]);}
#endif
#if (HIGH_PREC_TYPE==GMP_HIGH_PREC)
  int mpf_precision;
  int high_prec_nbits(){return mpf_get_default_prec();}
  void float_high_prec_t_print(float_high_prec_t a){master_printf("workaround high_prec %lg\n",a.get_d());}
#endif
  
  //takes integer power
  float_high_prec_t float_high_prec_t_pow_int(float_high_prec_t in,int d)
  {
    float_high_prec_t out;
    
    //null case
    if(d==0) out=(float_high_prec_t)1.0;
    else
      //negative case
      if(d<0)
        {
          //compute inv and assign to out
          float_high_prec_t inv=out=1/in;
          //multiply the remaining numer of times
          for(int i=2;i<=-d;i++) out*=inv;
        }
    //positive case
      else
        {
          out=in;
          for(int i=2;i<=d;i++) out*=in;
	}
    
    return out;
  }
  
  //frac power
  float_high_prec_t float_high_prec_t_pow_int_frac(float_high_prec_t ext_in,int n,int d)
  {
    //take a copy
    float_high_prec_t in=ext_in;
    
    //compute by solving out^d=in^n=ref
    float_high_prec_t ref=float_high_prec_t_pow_int(in,n);
    
    //let's start from a reasonable approx
    float_high_prec_t out=pow(in.get_d(),(double)n/d);
    
    //(out+err)^d=in^n -> err=out*rel_err, rel_err=(ref/out^d-1)/d
    [[maybe_unused]] int iter=0;
    float_high_prec_t rel_residue;
    double tol=32*pow(2.0,-high_prec_nbits());
    do
      {
        //compute out^d
        float_high_prec_t outd=float_high_prec_t_pow_int(out,d);
	
        //compute relative error
        rel_residue=(ref/outd-1);
        
        //total err
        float_high_prec_t err=rel_residue*out/d;
        out+=err;
        //verbosity_lv3_master_printf("Iter %d rel_residue: %lg, tol: %lg\n",iter,fabs(rel_residue.get_d()),tol);
	
        iter++;
      }
    while(fabs(rel_residue.get_d())>tol);
    
    return out;
  }
}
