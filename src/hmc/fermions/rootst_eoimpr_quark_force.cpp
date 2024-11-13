#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/bench.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "hmc/backfield.hpp"
#include "inverters/staggered/cgm_invert_stD2ee_m2.hpp"
#include "new_types/su3.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //Compute the fermionic force the rooted staggered eoprec improved theory.
  //Of the result still need to be taken the TA and product with U
  //The approximation need to be already scaled, and must contain physical mass term
  void summ_the_rootst_eoimpr_quark_force(eo_ptr<quad_su3> F,eo_ptr<quad_su3> eo_conf,color* pf,eo_ptr<quad_u1> u1b,rat_approx_t* appr,double residue)
  {
    crash("reimplement");
    
    // const int nterms=appr->degree();
    
    // START_TIMING(quark_force_over_time,nquark_force_over);
    
    // //allocate each terms of the expansion
    // color **v_o=nissa_malloc("v_o",nterms,color*),**chi_e=nissa_malloc("chi_e",nterms,color*);
    // for(int iterm=0;iterm<nterms;iterm++)
    //   {
    // 	v_o[iterm]=nissa_malloc("v_o",locVolh+bord_volh,color);
    // 	chi_e[iterm]=nissa_malloc("chi_e",locVolh+bord_volh,color);
    //   }
    
    // //add the background fields
    // add_backfield_with_stagphases_to_conf(eo_conf,u1b);
    
    // //invert the various terms
    // STOP_TIMING(quark_force_over_time);
    // inv_stD2ee_m2_cgm_run_hm_up_to_comm_prec(chi_e,eo_conf,appr->poles.data(),nterms,1000000,residue,pf);
    // UNPAUSE_TIMING(quark_force_over_time);
    
    // ////////////////////
    
    // //summ all the terms performing appropriate elaboration
    // //possible improvement by communicating more borders together
    // for(int iterm=0;iterm<nterms;iterm++) apply_stDoe(v_o[iterm],eo_conf,chi_e[iterm]);
    
    // //remove the background fields
    // rem_backfield_with_stagphases_from_conf(eo_conf,u1b);
    
    // //communicate borders of v_o (could be improved...)
    // for(int iterm=0;iterm<nterms;iterm++) communicate_od_color_borders(v_o[iterm]);
    
    // //conclude the calculation of the fermionic force
    // for(int iterm=0;iterm<nterms;iterm++)
    //   {
    // 	const double weight=appr->weights[iterm];
	
    // 	NISSA_PARALLEL_LOOP(ieo,0,locVolh)
    // 	for(int mu=0;mu<NDIM;mu++)
    // 	  for(int ic1=0;ic1<NCOL;ic1++)
    // 	    for(int ic2=0;ic2<NCOL;ic2++)
    // 	      {
    // 		complex temp1,temp2;
		
    // 		//this is for ieo=EVN
    // 		unsafe_complex_conj2_prod(temp1,v_o[iterm][loceo_neighup[EVN][ieo][mu]][ic1],chi_e[iterm][ieo][ic2]);
    // 		unsafe_complex_prod(temp2,temp1,u1b[EVN][ieo][mu]);
    // 		complex_summ_the_prod_double(F[EVN][ieo][mu][ic1][ic2],temp2,weight*get_stagphase_of_lx(loclx_of_loceo[EVN][ieo],mu));
		
    // 		//this is for ieo=ODD
    // 		unsafe_complex_conj2_prod(temp1,chi_e[iterm][loceo_neighup[ODD][ieo][mu]][ic1],v_o[iterm][ieo][ic2]);
    // 		unsafe_complex_prod(temp2,temp1,u1b[ODD][ieo][mu]);
    // 		complex_subt_the_prod_double(F[ODD][ieo][mu][ic1][ic2],temp2,weight*get_stagphase_of_lx(loclx_of_loceo[ODD][ieo],mu));
    // 	      }
    // 	NISSA_PARALLEL_LOOP_END;
    //   }
    
    // //free
    // for(int iterm=0;iterm<nterms;iterm++)
    //   {
    // 	nissa_free(v_o[iterm]);
    // 	nissa_free(chi_e[iterm]);
    //   }
    
    // nissa_free(chi_e);
    // nissa_free(v_o);
    
    // STOP_TIMING(quark_force_over_time);
  }
}
