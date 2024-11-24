/////////////////////////// MATRIX ELEMENT METHOD 4 ELLESETTETE ////////////////////////////
// We need to compute O=<0|P^0|pi_0> expanding the VEV around the Isosymmetric point up to 1st order in Δm
// Thus <O> = <O>_0 - ΔmΣ_x<0|P^0 L_ib(x)|pi_0>_0
// where ψ=(u,d),	L_ib(x)=ψ̅ σ_3ψ(x),	  P^0 = (u̅ g5 u + d̅ g5 d), 	pi_0 = (u̅ g5 d - d̅ g5 u)  *(Neglecting norm fact.)*
// It turns out that only diagrams that contribute to <O> are the following ones:
// (just considering one flavour)
/** 
	1.a)connected with insertion:
	   sink o(m)------->source o'(0)	   X:insertion of L_ib(n)    G(_|_)= D^-1(_|_) : porpagator isosymmetric point
													  
                      ***\X/ ***				Σ_n g5 g5 <0| u̅(m)u(m) u̅(n)u(n) u̅(0)u(0) |0> =
		   **           **				Σ_n g5 g5 <0| G(m|n) G(n|0) G(0|m) |0> = 				       		  		
		  o               o'	       ===>		Σ_n Tr[<0| G(m|n) G(n|0) g5 G(0|m) g5 |0>] =
		   **           **				Σ_n Tr[<0| G(m|n) G(n|0) G(m|0)^† |0>] = (TF and proj to 0 momentum)
		      *********					~Σ_n Σ_m Tr[<0| G(m|n) G(n|0) G(m|0)^†|0>]	

	1.b)connected without insertion(isoSym):
													  
		      *********				         g5 g5 <0| u̅(m)u(m) u̅(0)u(0) |0> =
		   **           **				 g5 g5 <0| G(m|0) G(0|m)  |0> = 				       		  		
		  o               o'	       ===>		 Tr[<0| G(m|0) G(m|0)^† |0>] = (TF and proj to 0 momentum)
		   **           **				~ Σ_m Tr[<0|G(m|0) G(m|0)^† |0>] 
		      *********

	2)disconnected:
		  ******        ******				Σ_n g5 g5 <0| u̅(m)u(m) u̅(n)u(n) u̅(0)u(0) |0> =
		*        *    *        *			Σ_n g5 g5 <0| G(m|n) G(n|m) G(0|0) |0> =
	       o          X  *          o'     ===>             Σ_n Tr[<0| G(m|n) G(n|m) g5 |0>] Tr[<0| G(0|0) g5 |0>] = (TF and proj to 0 momentum)
		*        *    *        *			~Σ_n Σ_m Tr[<0|G(m|n) G(n|m) g5|0>] Tr[<0| G(0|0)g5 |0>]
		  ******        ******				

*/


// Method (just a sketch) for the connected part with one-end trick:
/* 
	General idea:  given a naive 2pt function ~ Tr <  G(m|n) G(m|n)^† >  and an ensemble of Nr stochastic sources {η} s.t. 1/N Σ η η* = id...
	then one introduces a "φ" prop. φr(m) = Σn G(m|n) ηr(n)  which is sol to Σn D(m|n) φr(n) = ηr(m)
	turns out that (1/Nr) Σr φr(m) ηr(n)* is unbiased estimator of G(m|n) and the prod of twofem gives an estimate for the whole diagram

	BUTT we have an insertion i.e. Σn Tr[<0| G(m|n) G(n|0) G(m|0)^† |0>] 
	then we introduce a sequential "Φ" prop. Φr(m) = Σn S(m|n) ηr(n)  where S(m|n) is sol to Σn D(m|n) S(n|0) ~ G(m|0) id..
	thus Φ is sol to Σn D(m|n) Φ(n) = φ(m)
	So the diagram should be prod of one "normal" φ and the sequential Φ*

	one-end trick autmatically project to 0 momentum 
	
*/


#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "new_types/su3.hpp"

#include "stag.hpp"
#include "ellesettete.hpp"

namespace nissa
{
  using namespace stag;
  
  // correlators for the measurement of ellesettete
  void measure_ellesettete(eo_ptr<quad_su3> conf,theory_pars_t &theory_pars,ellesettete_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    crash("reimplement");
    // int nflavs=theory_pars.nflavs();
    
    // //open the file, allocate point result and source
    // FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    // complex *point_result=nissa_malloc("point_result",locVol,complex);
    // NEW_FIELD_T(source);
    // NEW_FIELD_T(g5_id_source);
    // NEW_FIELD_T(id_g5_source);
    
    // //vectors for propagators calculation
    // NEW_FIELD_T(SIMPLE_PROP);
    // NEW_FIELD_T(PROP_ID_G5);
    // NEW_FIELD_T(ID_G5_PROP_ID_G5);
    // NEW_FIELD_T(SEQ_PROP);
    // NEW_FIELD_T(SEQ_PROP_ID_G5);
    // NEW_FIELD_T(ID_G5_SEQ_PROP_ID_G5);
	
    
    // for(int icopy=0;icopy<meas_pars.ncopies;icopy++)
    //   {
    // 	//print conf id and copy id
    // 	master_fprintf(file," # %d\t%d\t\n",iconf,icopy);
	
    // 	//loop over flavors
    // 	for(int iflav=0;iflav<nflavs;iflav++)
    // 	  {
    // 	    if(theory_pars.quarks[iflav].discretiz!=ferm_discretiz::ROOT_STAG) crash("not defined for non-staggered quarks");
	    
    // 	    for(int glb_t=0;glb_t<glbSize[0];glb_t++)
    // 	      {
    // 		//vectors for output
    // 		NEW_TRACE_RES_VEC(Tr_two_pts_iso,glbSize[0]);
    // 		NEW_TRACE_RES_VEC(Tr_two_pts,glbSize[0]);
    // 		NEW_TRACE_RES_VEC(Tr_three_pts,glbSize[0]);
    // 		NEW_TRACE_RES_VEC(Tr_four_pts,glbSize[0]);
    // 		NEW_TRACE_RES_VEC(Tr_insertion_bubble,glbSize[0]);
    // 		NEW_TRACE_RES_VEC(Tr_no_insertion_bubble,glbSize[0]);
		
    // 		master_fprintf(file," # source time %d\n",glb_t);
		
    // 		//loop over hits
    // 		for(int ihit=0;ihit<meas_pars.nhits;ihit++)
    // 		  {
    // 		    //prepare the sources with the right structure spin x taste
    // 		    fill_source(source,glb_t,meas_pars.rnd_type);
    // 		    apply_stag_op(g5_id_source,conf,theory_pars.backfield[iflav],GAMMA_INT::GAMMA_5,GAMMA_INT::IDENTITY,source);
    // 		    apply_stag_op(id_g5_source,conf,theory_pars.backfield[iflav],GAMMA_INT::IDENTITY,GAMMA_INT::GAMMA_5,source);

    // 		    //compute std 2pts propagator G(m|n) ~ [D^-1(m|y) source(y)] source(n)* and simple sequential propagator
    // 		    MINV(SIMPLE_PROP,iflav,source);
    // 		    MINV(SEQ_PROP,iflav,SIMPLE_PROP);

    // 		    //compute  2pts propagator with id x g5 at source and apply id x g5 at sink
    // 		    MINV(PROP_ID_G5,iflav,id_g5_source);
    // 		    apply_stag_op(ID_G5_PROP_ID_G5,conf,theory_pars.backfield[iflav],GAMMA_INT::IDENTITY,GAMMA_INT::GAMMA_5,PROP_ID_G5);

    // 		    //compute sequential propagator with id x g5 at source and apply id x g5 at sink
    // 		    MINV(SEQ_PROP_ID_G5,iflav,PROP_ID_G5);
    // 		    apply_stag_op(ID_G5_SEQ_PROP_ID_G5,conf,theory_pars.backfield[iflav],GAMMA_INT::IDENTITY,GAMMA_INT::GAMMA_5,SEQ_PROP_ID_G5);
			

    // 		    //then glb reduction to compute the trace for the connected 2pts_iso, 2pts, 3pts and 4pts diagrams
    // 		    SUMM_THE_TIME_TRACE_PRINT_AT_LAST_HIT(Tr_two_pts_iso,SIMPLE_PROP,SIMPLE_PROP);
    // 		    SUMM_THE_TIME_TRACE_PRINT_AT_LAST_HIT(Tr_two_pts,SIMPLE_PROP,ID_G5_PROP_ID_G5);
    // 		    SUMM_THE_TIME_TRACE_PRINT_AT_LAST_HIT(Tr_three_pts,SIMPLE_PROP,ID_G5_SEQ_PROP_ID_G5);
    // 		    SUMM_THE_TIME_TRACE_PRINT_AT_LAST_HIT(Tr_four_pts,SEQ_PROP,ID_G5_SEQ_PROP_ID_G5);

    // 		    //////// disconnected ////////
    // 		    //here we need just simple seq prop with nothing at source
    // 		    if(ihit==meas_pars.nhits-1) master_fprintf(file," # Tr_no_insertion_bubble source time %d\n", glb_t);
    // 		    SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_no_insertion_bubble,g5_id_source,SIMPLE_PROP);
    // 		    if(ihit==meas_pars.nhits-1) master_fprintf(file,"\n # Tr_insertion_bubble source time %d\n", glb_t);
    // 		    SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_insertion_bubble,g5_id_source,SEQ_PROP);
    // 		    master_fprintf(file,"\n");
    // 		  }
    // 		  master_fprintf(file,"\n");
    // 	      }
    // 	  }
	
    // 	master_fprintf(file,"\n");
    //   }
    
    // //deallocate and close file
    // DELETE_FIELD_T(SIMPLE_PROP);
    // DELETE_FIELD_T(PROP_ID_G5);
    // DELETE_FIELD_T(ID_G5_PROP_ID_G5);
    // DELETE_FIELD_T(SEQ_PROP);
    // DELETE_FIELD_T(SEQ_PROP_ID_G5);
    // DELETE_FIELD_T(ID_G5_SEQ_PROP_ID_G5);
    // DELETE_FIELD_T(source);
    // DELETE_FIELD_T(g5_id_source);
    // DELETE_FIELD_T(id_g5_source);
    // nissa_free(point_result);
    // close_file(file);
  }
  
  //print
  std::string ellesettete_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasElleSettete\n";
    os<<base_fermionic_meas_t::get_str(full);
    
    return os.str();
  }
}
