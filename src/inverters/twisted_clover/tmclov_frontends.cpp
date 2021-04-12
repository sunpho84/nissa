#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "cgm_invert_tmclovQ2.hpp"

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/tmclovQ/reconstruct_tmclov_doublet.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"

namespace nissa
{
  //invert a set of propagators using the passed source
  //the output is stored in twisted basis, assuming that prop=su3spinspin[2][nmass][>=loc_vol]
  void compute_su3spinspin_tmclov_propagators_multi_mass(su3spinspin ***prop,quad_su3 *conf,double kappa,clover_term_t *Cl,double *mass,int nmass,int niter_max,double *req_res,su3spinspin *source)
  {
    //allocate temporary source
    spincolor *temp_source=nissa_malloc("temp_source",locVolWithBord.nastyConvert(),spincolor);
    //allocate temp_vec
    spincolor *temp_vec[2];
    temp_vec[0]=nissa_malloc("temp_vec[0]",locVolWithBord.nastyConvert(),spincolor);
    temp_vec[1]=nissa_malloc("temp_vec[1]",locVolWithBord.nastyConvert(),spincolor);
    //allocate nmass spincolors, for the cgm solutions
    spincolor **cgm_solution;
    cgm_solution=nissa_malloc("cgm_solution",nmass,spincolor*);
    for(int imass=0;imass<nmass;imass++) cgm_solution[imass]=nissa_malloc("cgm_solution[imass]",locVolWithBord.nastyConvert(),spincolor);
    
    //loop over the source dirac and color index
    for(int id=0;id<NDIRAC;id++)
      for(int ic=0;ic<NCOL;ic++)
	{
	  NISSA_LOC_VOL_LOOP(ivol)
	    get_spincolor_from_su3spinspin(temp_source[ivol.nastyConvert()],source[ivol.nastyConvert()],id,ic);
	  set_borders_invalid(temp_source);
	  
	  double init_time=take_time();
	  inv_tmclovDQ_cgm(cgm_solution,conf,kappa,Cl,mass,nmass,niter_max,req_res,temp_source);
	  master_printf("Finished the inversion of D*Q, dirac index %d, color %d in %g sec\n",id,ic,take_time()-init_time);
	  
	  //reconstruct the doublet
	  for(int imass=0;imass<nmass;imass++)
	    {
	      reconstruct_tmclov_doublet(temp_vec[0],temp_vec[1],conf,kappa,Cl,mass[imass],cgm_solution[imass]);
	      
	      //convert the id-th spincolor into the su3spinspin
	      for(int r=0;r<2;r++)
		{
		  NISSA_LOC_VOL_LOOP(ivol)
		    put_spincolor_into_su3spinspin(prop[r][imass][ivol.nastyConvert()],temp_vec[r][ivol.nastyConvert()],id,ic);
		  set_borders_invalid(prop[r]);
		}
	      verbosity_lv2_master_printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
	    }
	}
    
    //free temp vec
    nissa_free(temp_vec[1]);nissa_free(temp_vec[0]);
    nissa_free(temp_source);
    for(int imass=0;imass<nmass;imass++) nissa_free(cgm_solution[imass]);
    nissa_free(cgm_solution);
  }
}
