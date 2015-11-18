#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/random.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "geometry/geometry_eo.hpp"
#include "hmc/backfield.hpp"
#include "inverters/staggered/cg_invert_stD.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "routines/math_routines.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //multiply by M^-1
  THREADABLE_FUNCTION_7ARG(mult_Minv, color**,prop, quad_su3**,conf, quad_u1**,u1b, double,m, double,residue, color**,source, bool,phase_present)
  {
    if(!phase_present) addrem_stagphases_to_eo_conf(conf);
    add_backfield_to_conf(conf,u1b);
    inv_stD_cg(prop,conf,m,100000,residue,source);
    rem_backfield_from_conf(conf,u1b);
    if(!phase_present) addrem_stagphases_to_eo_conf(conf);
  }
  THREADABLE_FUNCTION_END
  void mult_Minv(color **prop,quad_su3 **conf,theory_pars_t *pars,int iflav,double residue,color **source,bool phase_present=false)
  {mult_Minv(prop,conf,pars->backfield[iflav],pars->quark_content[iflav].mass,residue,source,phase_present);}
  
  //compute the matrix element of the derivative of the dirac operator between two vectors
  //forward and backward derivative are stored separately, for a reason
  void compute_fw_bw_der_mel(complex *res_fw_bw,color **left,quad_su3 **conf,int mu,color **right,complex *point_result)
  {
    GET_THREAD_ID();
    
    color **right_fw_bw[2]={right,left};
    
    for(int fw_bw=0;fw_bw<2;fw_bw++)
      {
	vector_reset(point_result);
	for(int par=0;par<2;par++)
	  NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	    {
	      color v;
	      unsafe_su3_prod_color(v,conf[par][ieo][mu],right_fw_bw[fw_bw][!par][loceo_neighup[par][ieo][mu]]);
	      complex t;
	      if(fw_bw==0) color_scalar_prod(t,right_fw_bw[!fw_bw][par][ieo],v);
	      else         color_scalar_prod(t,v,right_fw_bw[!fw_bw][par][ieo]);
	      complex_summassign(point_result[loclx_of_loceo[par][ieo]],t);
	    }
	THREAD_BARRIER();
	complex_vector_glb_collapse(res_fw_bw[fw_bw],point_result,loc_vol);
      }
  }
}
