#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/bench.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec.hpp"
#include "dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec.hpp"
#include "hmc/backfield.hpp"
#include "inverters/twisted_clover/cgm_invert_tmclovDkern_eoprec_square_portable.hpp"
#include "new_types/su3.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  THREADABLE_FUNCTION_9ARG(summ_the_roottm_clov_eoimpr_quark_force, quad_su3**,F, quad_su3**,eo_conf, double,kappa, clover_term_t*,Cl_odd, inv_clover_term_t*,invCl_evn, spincolor*,phi_o, quad_u1**,u1b, rat_approx_t*,appr, double,residue)
  {
    GET_THREAD_ID();
    
    START_TIMING(quark_force_over_time,nquark_force_over);
    
    //allocate each terms of the expansion
    spincolor *Y[2][appr->degree()],*X[2][appr->degree()],*temp=nissa_malloc("temp",loc_volh+bord_volh,spincolor);
    for(int eo=0;eo<2;eo++)
      for(int iterm=0;iterm<appr->degree();iterm++)
	{
	  Y[eo][iterm]=nissa_malloc("Y",loc_volh+bord_volh,spincolor);
	  X[eo][iterm]=nissa_malloc("X",loc_volh+bord_volh,spincolor);
      }
    
    //add the background fields
    add_backfield_without_stagphases_to_conf(eo_conf,u1b);
    
    //invert the various terms
    STOP_TIMING(quark_force_over_time);
    inv_tmclovDkern_eoprec_square_portable_run_hm_up_to_comm_prec(X[ODD],eo_conf,kappa,Cl_odd,invCl_evn,appr->poles.data(),appr->degree(),10000000,residue,phi_o);
    UNPAUSE_TIMING(quark_force_over_time);
    
    ////////////////////
    
    // now we need to implement eq. B.10 of https://arxiv.org/pdf/0905.3331.pdf
    
    //summ all the terms performing appropriate elaboration
    //possible improvement by communicating more borders together
    for(int iterm=0;iterm<appr->degree();iterm++)
      {
	tmclovDkern_eoprec_eos(Y[ODD][iterm],temp,eo_conf,kappa,Cl_odd,invCl_evn,true,appr->poles[iterm],X[ODD][iterm]);
	
	spincolor *v_o[]={X[ODD][iterm],Y[ODD][iterm]};
	spincolor *v_e[]={X[EVN][iterm],Y[EVN][iterm]};
	bool d[]={true,false};
	
	for(int i=0;i<2;i++)
	  {
	    tmn2Deo_eos(temp,eo_conf,v_o[i]);
	    inv_tmclovDee_or_oo_eos(v_e[i],invCl_evn,d[i],temp);
	    double_vector_prodassign_double((double*)(v_e[i]),0.5,loc_volh*sizeof(spincolor)/sizeof(double));
	  }
      }
    //remove the background fields
    rem_backfield_without_stagphases_from_conf(eo_conf,u1b);
    
    //communicate borders (could be improved...)
    for(int eo=0;eo<2;eo++)
      for(int iterm=0;iterm<appr->degree();iterm++)
	for(auto v : {X[eo][iterm],Y[eo][iterm]})
	  communicate_ev_or_od_spincolor_borders(v,eo);
    
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
    for(int eo=0;eo<2;eo++)
      for(int iterm=0;iterm<appr->degree();iterm++)
	{
	  nissa_free(X[eo][iterm]);
	  nissa_free(Y[eo][iterm]);
	}
    nissa_free(temp);
    
    STOP_TIMING(quark_force_over_time);
  }
  THREADABLE_FUNCTION_END
}
