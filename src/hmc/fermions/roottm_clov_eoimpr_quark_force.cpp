#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/bench.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "hmc/backfield.hpp"
#include "inverters/twisted_clover/cgm_invert_tmclovDkern_eoprec_square_portable.hpp"
#include "new_types/su3.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  THREADABLE_FUNCTION_10ARG(summ_the_roottm_clov_eoimpr_quark_force, quad_su3**,F, quad_su3**,eo_conf, double,kappa, clover_term_t*,Cl_odd, inv_clover_term_t*,invCl_evn, double,mass, spincolor*,pf, quad_u1**,u1b, rat_approx_t*,appr, double,residue)
  {
    GET_THREAD_ID();
    
    START_TIMING(quark_force_over_time,nquark_force_over);
    
    //allocate each terms of the expansion
    spincolor *v_o[appr->degree()],*chi_e[appr->degree()];
    for(int iterm=0;iterm<appr->degree();iterm++)
      {
	v_o[iterm]=nissa_malloc("v_o",loc_volh+bord_volh,spincolor);
	chi_e[iterm]=nissa_malloc("chi_e",loc_volh+bord_volh,spincolor);
      }
    
    //add the background fields
    add_backfield_without_stagphases_to_conf(eo_conf,u1b);
    
    //invert the various terms
    STOP_TIMING(quark_force_over_time);
    inv_tmclovDkern_eoprec_square_portable_run_hm_up_to_comm_prec(chi_e,eo_conf,kappa,Cl_odd,invCl_evn,mass,appr->poles.data(),appr->degree(),10000000,residue,pf);
    //inv_stD2ee_m2_cgm_run_hm_up_to_comm_prec(chi_e,eo_conf,appr->poles.data(),appr->degree(),1000000,residue,pf);
    UNPAUSE_TIMING(quark_force_over_time);
    
    ////////////////////
    
    //summ all the terms performing appropriate elaboration
    //possible improvement by communicating more borders together
    crash("");
    //for(int iterm=0;iterm<appr->degree();iterm++) apply_stDoe(v_o[iterm],eo_conf,chi_e[iterm]);
    
    //remove the background fields
    rem_backfield_without_stagphases_from_conf(eo_conf,u1b);
    
    //communicate borders of v_o (could be improved...)
    crash("");
    //for(int iterm=0;iterm<appr->degree();iterm++) communicate_od_color_borders(v_o[iterm]);
    
    //conclude the calculation of the fermionic force
    for(int iterm=0;iterm<appr->degree();iterm++)
      NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	for(int mu=0;mu<NDIM;mu++)
	  for(int ic1=0;ic1<NCOL;ic1++)
	    for(int ic2=0;ic2<NCOL;ic2++)
	      {
		crash("");
		// complex temp1,temp2;
		
		// //this is for ieo=EVN
		// unsafe_complex_conj2_prod(temp1,v_o[iterm][loceo_neighup[EVN][ieo][mu]][ic1],chi_e[iterm][ieo][ic2]);
		// unsafe_complex_prod(temp2,temp1,u1b[EVN][ieo][mu]);
		// complex_summ_the_prod_double(F[EVN][ieo][mu][ic1][ic2],temp2,appr->weights[iterm]*get_stagphase_of_lx(loclx_of_loceo[EVN][ieo],mu));
		
		// //this is for ieo=ODD
		// unsafe_complex_conj2_prod(temp1,chi_e[iterm][loceo_neighup[ODD][ieo][mu]][ic1],v_o[iterm][ieo][ic2]);
		// unsafe_complex_prod(temp2,temp1,u1b[ODD][ieo][mu]);
		// complex_subt_the_prod_double(F[ODD][ieo][mu][ic1][ic2],temp2,appr->weights[iterm]*get_stagphase_of_lx(loclx_of_loceo[ODD][ieo],mu));
	      }
    NISSA_PARALLEL_LOOP_END;
    
    //free
    for(int iterm=0;iterm<appr->degree();iterm++)
      {
	nissa_free(v_o[iterm]);
	nissa_free(chi_e[iterm]);
      }
    
    STOP_TIMING(quark_force_over_time);
  }
  THREADABLE_FUNCTION_END
}
