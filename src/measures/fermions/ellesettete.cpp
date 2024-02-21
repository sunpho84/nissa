/////////////////////////// MATRIX ELEMENT METHOD 4 ELLESETTETE ////////////////////////////
// We need to compute O=<0|P^0|pi_0> expanding the VEV around the Isosymmetric point up to 1st order in Δm
// Thus <O> = <O>_0 - ΔmΣ_x<0|P^0 L_ib(x)|pi_0>_0
// where ψ=(u,d),	L_ib(x)=ψ̅ σ_3ψ(x),	  P^0 = (u̅ g5 u + d̅ g5 d), 	pi_0 = (u̅ g5 d - d̅ g5 u)  *(Neglecting norm fact.)*
// It turns out that only diagrams that contribute to <O> are the following ones:
// (just considering one flavour)
/* 
	1)connected:
	   sink o(m)------->source o'(0)	   X:insertion of L_ib(n)    G(_|_)= D^-1(_|_) : porpagator isosymmetric point
													  
			  

                      ***\X/***					Σ_n g5 g5 <0| u̅(m)u(m) u̅(n)u(n) u̅(0)u(0) |0> =
		   **           **				Σ_n g5 g5 <0| G(m|n) G(n|0) G(0|m) |0> = 				       		  		
		  o               o'	       ===>		Σ_n Tr[<0| G(m|n) G(n|0) g5 G(0|m) g5 |0>] =
		   **           **				Σ_n Tr[<0| G(m|n) G(n|0) G(m|0)^† |0>] = (TF and proj to 0 momentum)
		      *********					~Σ_n Σ_m Tr[<0| G(m|n) G(n|0) G(m|0)^†|0>]	


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

	one-end trick autmatically project to 0 momentum (?)

*/


#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"

#include "stag.hpp"
#include "ellesettete.hpp"

namespace nissa
{
  using namespace stag;
  
  // measure the chiral condensate and its derivative w.r.t mu
  void measure_ellesettete(eo_ptr<quad_su3> conf,theory_pars_t &theory_pars,ellesettete_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    int nflavs=theory_pars.nflavs();
    
    //open the file, allocate point result and source
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    complex *point_result=nissa_malloc("point_result",locVol,complex);
    NEW_FIELD_T(source);
    NEW_FIELD_T(g5_source);//----> not needed.

    //vectors for calculation
    NEW_FIELD_T(SIMPLE_PROP);
    NEW_FIELD_T(SEQ_PROP);
    //NEW_FIELD_T(G5_PROP);
   
    
    for(int icopy=0;icopy<meas_pars.ncopies;icopy++)
      {
	//print conf id and copy id
	master_fprintf(file,"%d\t%d\t",iconf,icopy);
	
	//loop over flavors
	for(int iflav=0;iflav<nflavs;iflav++)
	  {
	    if(theory_pars.quarks[iflav].discretiz!=ferm_discretiz::ROOT_STAG) crash("not defined for non-staggered quarks");
	    for(int glb_t=0; glb_t<glbSize[0]; glb_t++)
		{
			//vectors for output
			NEW_TRACE_RES(Tr_three_pts);
			NEW_TRACE_RES_VEC(Tr_first_bubble, glbSize[0]);
			NEW_TRACE_RES_VEC(Tr_second_bubble, glbSize[0]);
						
	    
	    	//loop over hits
			for(int ihit=0;ihit<meas_pars.nhits;ihit++)
			{
				//fill the source
				fill_source(source,glb_t,meas_pars.rnd_type);
				
				//compute std 2pts propagator G(m|n) ~ [D^-1(m|y) source(y)] source(n)*
				MINV(SIMPLE_PROP,iflav,source);
				
				//compute sequential propagator G(m|n) ~ [D^-1(m|y) source(y)] D^-1(y|n)
				MINV(SEQ_PROP,iflav,SIMPLE_PROP);

				//then glb reduction to compute the trace for the connected 3pts diagram 
				SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_three_pts,SIMPLE_PROP,SEQ_PROP);
				
				//////// disconnected //////// 
				apply_stag_op(g5_source,conf,theory_pars->backfield[iflav],source);
				
				//MINV(G5_PROP,iflav,g5_source); 
				SUMM_THE_TIME_TRACE_PRINT_AT_LAST_HIT(Tr_first_bubble, SEQ_PROP, g5_source);
				SUMM_THE_TIME_TRACE_PRINT_AT_LAST_HIT(Tr_second_bubble,SIMPLE_PROP,15,15,g5_source);
			}
		}	
	  }
	
	master_fprintf(file,"\n");
      }
    
    //deallocate and close file
    DELETE_FIELD_T(SIMPLE_PROP);
    DELETE_FIELD_T(SEQ_PROP);
    DELETE_FIELD_T(G5_PROP);
    
    
    close_file(file);
    nissa_free(point_result);
    DELETE_FIELD_T(source);
	DELETE_FIELD_T(g5_source);
  }
  
  //print
  std::string ellesettete_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasElle7\n";
    os<<base_fermionic_meas_t::get_str(full);
    
    return os.str();
  }
}
