#ifndef _CLOVER_TERM_BGQ_HPP
#define _CLOVER_TERM_BGQ_HPP

#include "bgq/bgq_macros.hpp"
#include "bgq/intrinsic.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  inline void unsafe_vir_halfspincolor_halfspincolor_times_vir_halfspincolor(vir_halfspincolor out,vir_halfspincolor_halfspincolor invCl,vir_halfspincolor in)
  {
    DECLARE_REG_VIR_HALFSPINCOLOR(REG_IN);
    REG_LOAD_VIR_HALFSPINCOLOR(REG_IN,in);
    
    for(int id=0;id<NDIRAC/2;id++)
      for(int ic=0;ic<NCOL;ic++)
	{
	  DECLARE_REG_VIR_HALFSPINCOLOR(M);
	  REG_LOAD_VIR_HALFSPINCOLOR(M,invCl[id][ic]);
	  
	  DECLARE_REG_VIR_COMPLEX(REG_OUT);
	  
	  REG_VIR_COMPLEX_PROD(REG_OUT,M_s0_c0,REG_IN_s0_c0);
	  REG_VIR_COMPLEX_SUMM_THE_PROD(REG_OUT,M_s0_c1,REG_IN_s0_c1);
	  REG_VIR_COMPLEX_SUMM_THE_PROD(REG_OUT,M_s0_c2,REG_IN_s0_c2);
	  REG_VIR_COMPLEX_SUMM_THE_PROD(REG_OUT,M_s1_c0,REG_IN_s1_c0);
	  REG_VIR_COMPLEX_SUMM_THE_PROD(REG_OUT,M_s1_c1,REG_IN_s1_c1);
	  REG_VIR_COMPLEX_SUMM_THE_PROD(REG_OUT,M_s1_c2,REG_IN_s1_c2);
	  
	  STORE_REG_VIR_COMPLEX(out[id][ic],REG_OUT);
	}
  }
  
  inline void unsafe_vir_halfspincolor_halfspincolor_dag_times_vir_halfspincolor(vir_halfspincolor out,vir_halfspincolor_halfspincolor invCl,vir_halfspincolor in)
  {
    DECLARE_REG_VIR_HALFSPINCOLOR(REG_OUT);
    REG_SPLAT_VIR_HALFSPINCOLOR(REG_OUT,0);
    
    for(int id=0;id<NDIRAC/2;id++)
      for(int ic=0;ic<NCOL;ic++)
	{
	  DECLARE_REG_VIR_COMPLEX(REG_IN);
	  REG_LOAD_VIR_COMPLEX(REG_IN,in[id][ic]);
	  
	  DECLARE_REG_VIR_HALFSPINCOLOR(M);
	  REG_LOAD_VIR_HALFSPINCOLOR(M,invCl[id][ic]);
	  
	  REG_VIR_COMPLEX_SUMM_THE_CONJ1_PROD(REG_OUT_s0_c0,M_s0_c0,REG_IN);
	  REG_VIR_COMPLEX_SUMM_THE_CONJ1_PROD(REG_OUT_s0_c1,M_s0_c1,REG_IN);
	  REG_VIR_COMPLEX_SUMM_THE_CONJ1_PROD(REG_OUT_s0_c2,M_s0_c2,REG_IN);
	  REG_VIR_COMPLEX_SUMM_THE_CONJ1_PROD(REG_OUT_s1_c0,M_s1_c0,REG_IN);
	  REG_VIR_COMPLEX_SUMM_THE_CONJ1_PROD(REG_OUT_s1_c1,M_s1_c1,REG_IN);
	  REG_VIR_COMPLEX_SUMM_THE_CONJ1_PROD(REG_OUT_s1_c2,M_s1_c2,REG_IN);
	}
    
    STORE_REG_VIR_HALFSPINCOLOR(out,REG_OUT);
  }
}

#endif
