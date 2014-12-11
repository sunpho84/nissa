#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/macros.hpp"
#include "base/thread_macros.hpp"
#include "bgq/bgq_macros.hpp"
#include "bgq/Wilson_hopping_matrix_eo_or_oe_bgq.hpp"
//debug
#include "geometry/geometry_vir.hpp"
#include "routines/ios.hpp"

#include "new_types/new_types_definitions.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  void inv_tmDee_or_oo_eos(bi_spincolor *out,double kappa,double mass,bi_spincolor *in)
  {
    double a=1/(2*kappa),b=mass,nrm=a*a+b*b;
    a/=nrm;
    b/=-nrm;
    bi_complex z={a,b,a,b};
    DECLARE_REG_BI_COMPLEX(BI_REG_Z);
    REG_LOAD_BI_COMPLEX(BI_REG_Z,z);
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(X,0,loc_volh/2)
      {
	DECLARE_REG_BI_HALFSPINCOLOR(BI_REG_IN);
	DECLARE_REG_BI_HALFSPINCOLOR(BI_REG_OUT);
	
	REG_LOAD_BI_HALFSPINCOLOR(BI_REG_IN,in[X][0]);
	REG_BI_HALFSPINCOLOR_PROD_COMPLEX(BI_REG_OUT,BI_REG_IN,BI_REG_Z);
	STORE_REG_BI_HALFSPINCOLOR(out[X][0],BI_REG_OUT);
	
	REORDER_BARRIER();
	
	REG_LOAD_BI_HALFSPINCOLOR(BI_REG_IN,in[X][2]);
	REG_BI_HALFSPINCOLOR_PROD_CONJ2_COMPLEX(BI_REG_OUT,BI_REG_IN,BI_REG_Z);
	STORE_REG_BI_HALFSPINCOLOR(out[X][2],BI_REG_OUT);
      }
    
    set_borders_invalid(out);
  }
  
  //implement ee or oo part of Dirac operator, equation(3) and finish
  void minus_one_quarter_g5_subt_tmQee_or_oo_eos(bi_spincolor *out,double kappa,double mass,bi_spincolor *in)
  {
    DECLARE_REG_BI_COMPLEX(BI_REG_MONE_QUARTER);
    REG_SPLAT_BI_COMPLEX(BI_REG_MONE_QUARTER,-0.25);
    DECLARE_REG_BI_COMPLEX(BI_REG_ONE_QUARTER);
    REG_SPLAT_BI_COMPLEX(BI_REG_ONE_QUARTER,+0.25);
    
    double a=1/(2*kappa),b=mass;
    bi_complex z={a,b, a,b};
    bi_complex mzc={-a,b, -a,b};
    DECLARE_REG_BI_COMPLEX(BI_REG_Z);
    DECLARE_REG_BI_COMPLEX(BI_REG_MZC);
    REG_LOAD_BI_COMPLEX(BI_REG_Z,z);
    REG_LOAD_BI_COMPLEX(BI_REG_MZC,mzc);
    
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_volh/2)
      {
	DECLARE_REG_BI_HALFSPINCOLOR(BI_REG_OUT);
	DECLARE_REG_BI_HALFSPINCOLOR(BI_REG_IN);

	//multiply OUT{[0],[1]} by -1/4
	REG_LOAD_BI_HALFSPINCOLOR(BI_REG_OUT,out[ivol][0]);
	REG_BI_HALFSPINCOLOR_PROD_4DOUBLE(BI_REG_OUT,BI_REG_OUT,BI_REG_MONE_QUARTER);
	
	//summ Z*in
	REG_LOAD_BI_HALFSPINCOLOR(BI_REG_IN,in[ivol][0]);
	REG_BI_HALFSPINCOLOR_SUMM_THE_PROD_COMPLEX(BI_REG_OUT,BI_REG_IN,BI_REG_Z);
        STORE_REG_BI_HALFSPINCOLOR(out[ivol][0],BI_REG_OUT);
	
	REORDER_BARRIER();
	
	//multiply OUT{[2],[3]} by +1/4
	REG_LOAD_BI_HALFSPINCOLOR(BI_REG_OUT,out[ivol][2]);
	REG_BI_HALFSPINCOLOR_PROD_4DOUBLE(BI_REG_OUT,BI_REG_OUT,BI_REG_ONE_QUARTER);
	
	//subt (Z^+)*in
	REG_LOAD_BI_HALFSPINCOLOR(BI_REG_IN,in[ivol][2]);
	REG_BI_HALFSPINCOLOR_SUMM_THE_PROD_COMPLEX(BI_REG_OUT,BI_REG_IN,BI_REG_MZC);
        STORE_REG_BI_HALFSPINCOLOR(out[ivol][2],BI_REG_OUT);
      }
    set_borders_invalid(out);
  }
  
  void tmn2Doe_or_tmn2Deo_eos_bgq(bi_spincolor *out,bi_oct_su3 **conf,int oe_or_eo,bi_spincolor *in)
  {
    //compute on the surface and start communications
    apply_double_Wilson_hopping_matrix_oe_or_eo_bgq_nocomm(conf,0,vsurf_volh,in,oe_or_eo);
    start_double_Wilson_hopping_matrix_oe_or_eo_bgq_communications();
    
    //compute on the bulk and finish communications
    apply_double_Wilson_hopping_matrix_oe_or_eo_bgq_nocomm(conf,vsurf_volh,loc_volh/2,in,oe_or_eo);
    finish_double_Wilson_hopping_matrix_oe_or_eo_bgq_communications(oe_or_eo);
    
    //put together all the 8 pieces
    hopping_matrix_oe_or_eo_expand_to_double_Wilson_2D_bgq(out);
  }
  
  //wrappers
  void tmn2Doe_eos_bgq(bi_spincolor *out,bi_oct_su3 **conf,bi_spincolor *in){tmn2Doe_or_tmn2Deo_eos_bgq(out,conf,0,in);}
  void tmn2Deo_eos_bgq(bi_spincolor *out,bi_oct_su3 **conf,bi_spincolor *in){tmn2Doe_or_tmn2Deo_eos_bgq(out,conf,1,in);}
  
  //implement Koo defined in equation (7)
  void tmDkern_eoprec_eos_bgq(bi_spincolor *out,bi_oct_su3 **conf,double kappa,double mass,bi_spincolor *in)
  {
    if(in==out) crash("cannot work with in==out");
    
    tmn2Deo_eos_bgq(out,conf,in);
    inv_tmDee_or_oo_eos(out,kappa,mass,out);
    tmn2Doe_eos_bgq(out,conf,out);
    
    minus_one_quarter_g5_subt_tmQee_or_oo_eos(out,kappa,mass,in);
  }

  //square of Koo
  THREADABLE_FUNCTION_6ARG(tmDkern_eoprec_square_eos_bgq, bi_spincolor*,out, bi_spincolor*,temp, bi_oct_su3**,conf, double,kappa, double,mass, bi_spincolor*,in)
  {
    tmDkern_eoprec_eos_bgq(temp,conf,kappa,-mass, in);
    tmDkern_eoprec_eos_bgq(out,conf,kappa,+mass,temp);
  }
  THREADABLE_FUNCTION_END
}
