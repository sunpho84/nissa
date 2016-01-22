#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/complex.hpp"
#include "new_types/dirac.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3.hpp"

namespace nissa
{
  ///////////////////////// take Tr[g1 * s1^dag * g2 * s2], useful for mesons 2 points //////////////////
  
  void trace_g_ss_dag_g_ss(complex c,dirac_matr *g1,spinspin s1,dirac_matr *g2,spinspin s2)
  {
    spinspin t1,t2;
    
    unsafe_dirac_prod_spinspin_dag(t1,g1,s1);
    unsafe_dirac_prod_spinspin(t2,g2,s2);
    
    trace_prod_spinspins(c,t1,t2);
  }
  void trace_g_css_dag_g_css(complex c,dirac_matr *g1,colorspinspin s1,dirac_matr *g2,colorspinspin s2)
  {
    //reset out
    c[0]=c[1]=0;
    
    //loop over color indices
    for(int ic=0;ic<3;ic++)
      {
	spinspin t1,t2;
	
	unsafe_dirac_prod_spinspin_dag(t1,g1,s1[ic]);
	unsafe_dirac_prod_spinspin(t2,g2,s2[ic]);
	
	summ_the_trace_prod_spinspins(c,t1,t2);
      }
  }
  void trace_g_ccss_dag_g_ccss(complex c,dirac_matr *g1,su3spinspin s1,dirac_matr *g2,su3spinspin s2)
  {
    //reset out
    c[0]=c[1]=0;
    
    //loop over color indices
    for(int ic1=0;ic1<3;ic1++)
      for(int ic2=0;ic2<3;ic2++)
	{
	  spinspin t1,t2;
	  
	  unsafe_dirac_prod_spinspin_dag(t1,g1,s1[ic2][ic1]);
	  unsafe_dirac_prod_spinspin(t2,g2,s2[ic2][ic1]);
	  
	  summ_the_trace_prod_spinspins(c,t1,t2);
	}
  }
  
  ////////////////////////////////////////// BK contraction ///////////////////////////////////////
  
  void trace_g_ss_dag_g_ss_g_ss_dag_g_ss(complex c,dirac_matr *g1,spinspin s1,dirac_matr *g2,spinspin s2, dirac_matr *g3,spinspin s3,dirac_matr *g4,spinspin s4)
  {
    spinspin t1,t2,t12,t3,t4,t34;
    
    unsafe_dirac_prod_spinspin_dag(t1,g1,s1);
    unsafe_dirac_prod_spinspin(t2,g2,s2);
    unsafe_spinspin_prod_spinspin(t12,t1,t2);
    
    unsafe_dirac_prod_spinspin_dag(t3,g3,s3);
    unsafe_dirac_prod_spinspin(t4,g4,s4);
    unsafe_spinspin_prod_spinspin(t34,t3,t4);
    
    trace_prod_spinspins(c,t12,t34);
  }
}
