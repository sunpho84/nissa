#ifndef _DIRAC_OPERATOR_TMDEOIMPR_BGQ_HPP
#define _DIRAC_OPERATOR_TMDEOIMPR_BGQ_HPP

#include "bgq/Wilson_hopping_matrix_eo_or_oe_bgq.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  void tmD_or_Qkern_eoprec_eos_bgq(vir_spincolor *out,vir_oct_su3 **conf,double kappa,double mu,vir_spincolor *in,int D_or_Q);
  
  inline void tmn2Doe_or_tmn2Deo_eos_bgq(vir_spincolor *out,vir_oct_su3 **conf,int oe_or_eo,vir_spincolor *in)
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
  inline void tmn2Doe_eos_bgq(vir_spincolor *out,vir_oct_su3 **conf,vir_spincolor *in){tmn2Doe_or_tmn2Deo_eos_bgq(out,conf,0,in);}
  inline void tmn2Deo_eos_bgq(vir_spincolor *out,vir_oct_su3 **conf,vir_spincolor *in){tmn2Doe_or_tmn2Deo_eos_bgq(out,conf,1,in);}

  //non squared
  inline void tmDkern_eoprec_eos_bgq(vir_spincolor *out,vir_oct_su3 **conf,double kappa,double mu,vir_spincolor *in)
  {tmD_or_Qkern_eoprec_eos_bgq(out,conf,kappa,mu,in,0);}
  inline void tmQkern_eoprec_eos_bgq(vir_spincolor *out,vir_oct_su3 **conf,double kappa,double mu,vir_spincolor *in)
  {tmD_or_Qkern_eoprec_eos_bgq(out,conf,kappa,mu,in,1);}
  
  //squared
  inline void tmD_or_Qkern_eoprec_square_eos_bgq(vir_spincolor *out,vir_spincolor *temp,vir_oct_su3 **conf,double kappa,double mu,vir_spincolor *in,int D_or_Q)
  {
    tmD_or_Qkern_eoprec_eos_bgq(temp,conf,kappa,-mu, in   ,D_or_Q);
    tmD_or_Qkern_eoprec_eos_bgq(out,  conf,kappa,+mu, temp,D_or_Q);
  }
  inline void tmDkern_eoprec_square_eos_bgq(vir_spincolor *out,vir_spincolor *temp,vir_oct_su3 **conf,double kappa,double mu,vir_spincolor *in)
  {tmD_or_Qkern_eoprec_square_eos_bgq(out,temp,conf,kappa,mu,in,0);}
  inline void tmQkern_eoprec_square_eos_bgq(vir_spincolor *out,vir_spincolor *temp,vir_oct_su3 **conf,double kappa,double mu,vir_spincolor *in)
  {tmD_or_Qkern_eoprec_square_eos_bgq(out,temp,conf,kappa,mu,in,1);}
}

#endif
