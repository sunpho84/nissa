#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/su3.hpp"
#include "linalgs/linalgs.hpp"
#include "threads/threads.hpp"
#include "operations/su3_paths/clover_term.hpp"
#include "hmc/quark_pars.hpp"

#include <vector>

namespace nissa
{
  //Computes the clover determinant action
  THREADABLE_FUNCTION_3ARG(clover_det_action, double*,act, std::vector<quark_content_t>,quark_content, eo_ptr<quad_su3>,eo_conf)
  {
    double res=0.0;
    bool need=false;
	for(auto& q : quark_content)
	  need|=(q.cSW!=0);
    
    if(not need)
      res=0.0;
    else
      {
	GET_THREAD_ID();
	double *loc_act=nissa_malloc("loc_act",loc_volh,double);
	
	eo_ptr<clover_term_t> Cl={NULL,NULL};
	for(int eo=0;eo<2;eo++) Cl[eo]=nissa_malloc("Cl",loc_volh,clover_term_t);
	chromo_operator(Cl,eo_conf);
	
	for(auto& q : quark_content)
	  if(q.cSW)
	    {
	      chromo_operator_include_cSW(Cl,q.cSW);
	      
	      const double &mass=q.mass;
	      const double &kappa=q.kappa;
	      
	      NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
		{
		  complex d[2];
		  for(int x_high_low=0;x_high_low<2;x_high_low++)
		    {
		      halfspincolor_halfspincolor e;
		      
		      fill_point_twisted_clover_term(e,x_high_low,Cl[EVN][ieo],mass,kappa);
		      
		      matrix_determinant(d[x_high_low],(complex*)e,NDIRAC*NCOL/2);
		    }
		  
		  //Product of the two subblocks determinants
		  complex p;
		  unsafe_complex_prod(p,d[0],d[1]);
		  
		  loc_act[ieo]=log(complex_norm2(p));
		}
	      NISSA_PARALLEL_LOOP_END;
	      THREAD_BARRIER();
	      
	      chromo_operator_remove_cSW(Cl,q.cSW);
	      
	      double flav_act;
	      double_vector_glb_collapse(&flav_act,loc_act,loc_volh);
	      
	      //half volume, all colors, all dirac, norm2. Deg is included below
	      const double offset=log((1/sqr(2*q.kappa)+sqr(q.mass)))*NCOL*NDIRAC*glb_volh;
	      flav_act-=offset;
	      
	      res+=-flav_act*q.deg;
	    }
	
	nissa_free(loc_act);
	for(int eo=0;eo<2;eo++) nissa_free(Cl[eo]);
      }
    
    *act=res;
  }
  THREADABLE_FUNCTION_END;
}
