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
    
    //vectors for calculation
    NEW_FIELD_T(M);           // M^-1
    NEW_FIELD_T(dM_M);        // M' M^-1
    NEW_FIELD_T(d2M_M);       // M'' M^-1
    NEW_FIELD_T(M_M);         // M^-2
    NEW_FIELD_T(dM_M_M);      // M' M^-2
    NEW_FIELD_T(d2M_M_M);     // M'' M^-2
    NEW_FIELD_T(dM_M_dM_M);   // (M' M^-1)^2
    NEW_FIELD_T(M_dM_M_dM_M); // M^-1 (M' M^-1)^2
    NEW_FIELD_T(TMP);         // parking variable
    
    for(int icopy=0;icopy<meas_pars.ncopies;icopy++)
      {
	//print conf id and copy id
	master_fprintf(file,"%d\t%d\t",iconf,icopy);
	
	//loop over flavors
	for(int iflav=0;iflav<nflavs;iflav++)
	  {
	    if(theory_pars.quarks[iflav].discretiz!=ferm_discretiz::ROOT_STAG) crash("not defined for non-staggered quarks");
	    
	    //vectors for output
	    NEW_TRACE_RES(Tr_M);
	    NEW_TRACE_RES(Tr_dM_M);
	    NEW_TRACE_RES(Tr_d2M_M);
	    NEW_TRACE_RES(Tr_M_M);
	    NEW_TRACE_RES(Tr_dM_M_M);
	    NEW_TRACE_RES(Tr_d2M_M_M);
	    NEW_TRACE_RES(Tr_dM_M_dM_M);
	    NEW_TRACE_RES(Tr_M_dM_M_dM_M);
	    
	    //loop over hits
	    for(int ihit=0;ihit<meas_pars.nhits;ihit++)
	      {
		//fill the source
		fill_source(source,-1,meas_pars.rnd_type);
		
		//compute M^-1, M' M^-1, M'' M^-1
		MINV(M,iflav,source);
		DMDMU(dM_M,iflav,1,M);
		DMDMU(d2M_M,iflav,2,M);
		
		//compute M^-2, M' M^-2, M'' M^-2
		MINV(M_M,iflav,M);
		DMDMU(dM_M_M,iflav,1,M_M);
		DMDMU(d2M_M_M,iflav,2,M_M);
		
		//compute (M' M^-1)^2
		DMDMU(dM_M_dM_M,iflav,1,M); // M' M^-1
		MINV(TMP,iflav,dM_M_dM_M); // M^-1 M' M^-1
		DMDMU(dM_M_dM_M,iflav,1,TMP);
		
		//compute M^-1 (M' M^-1)^2
		MINV(M_dM_M_dM_M,iflav,dM_M_dM_M);
		
		//print traces
		SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_M,source,M);
		SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_dM_M,source,dM_M);
		SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_d2M_M,source,d2M_M);
		SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_M_M,source,M_M);
		SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_dM_M_M,source,dM_M_M);
		SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_d2M_M_M,source,d2M_M_M);
		SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_dM_M_dM_M,source,dM_M_dM_M);
		SUMM_THE_TRACE_PRINT_AT_LAST_HIT(Tr_M_dM_M_dM_M,source,M_dM_M_dM_M);
	      }
	  }
	
	master_fprintf(file,"\n");
      }
    
    //deallocate and close file
    DELETE_FIELD_T(M);
    DELETE_FIELD_T(dM_M);
    DELETE_FIELD_T(d2M_M);
    DELETE_FIELD_T(M_M);
    DELETE_FIELD_T(dM_M_M);
    DELETE_FIELD_T(d2M_M_M);
    DELETE_FIELD_T(dM_M_dM_M);
    DELETE_FIELD_T(M_dM_M_dM_M);
    DELETE_FIELD_T(TMP);
    
    close_file(file);
    nissa_free(point_result);
    DELETE_FIELD_T(source);
  }
  
void print_corr(eo_ptr<quad_su3> ext_conf,theory_pars_t &tp,ellesettete_meas_pars_t &meas_pars,int iconf,int conf_created){
	double norm=1.0/(meas_pars.nhits*glbSpatVol);
	int nflavs=tp.nflavs();

	int size_corr=glbSize[0]*nflavs;;

	complex *corr=nissa_malloc("corr", size_corr, complex);
	
	
	//maybe more generalizable? like for whatever spin/taste cobination for source/sink and insertion? boh 
	int ncopies=meas_pars.ncopies;  // ncopies is the number of gauge configurations??? idk yet (???), if not delete it (same below)
									// if it is, we can do just one loop over copies and compute both connected and disconnected at each it. 
	FILE *file1=open_file(path1,conf_created?"w":"a"); //i would use std::ofstream 
		for(int icopy=0;icopy<ncopies;icopy++){
			
			verbosity_lv2_master_printf("Computing copy %d/%d\n",icopy,ncopies);
			//compute_connected(corr,ext_conf,&tp,&meas_pars) -----> to be implemented 
			for(int iflav=0;iflav<nflavs;iflav++)
			{

				master_fprintf(file1," # conf %d ;"
						" flv = %d , m = %lg\n",
						iconf,iflav,tp.quarks[iflav].mass);

				for(int t=0;t<glbSize[0];t++)
				{
				int ic = iflav*glbSize[0]+t;
				master_fprintf(file1,"%d %+16.16lg %+16.16lg\n",t,corr[ic][RE]*norm,corr[ic][IM]*norm);
				}
				master_fprintf(file1,"\n");
			}

		}
	close_file(file1);

	FILE *file2=open_file(path2,conf_created?"w":"a");
		for(int icopy=0;icopy<ncopies;icopy++)
		{
			verbosity_lv2_master_printf("Computing copy %d/%d\n",icopy,ncopies);
			//compute_disconnected(corr,ext_conf,&tp,&meas_pars)-----> to be implemented 
			for(int iflav=0;iflav<nflavs;iflav++)
			{

				master_fprintf(file2," # conf %d ;"
						" flv = %d , m = %lg\n",
						iconf,iflav,tp.quarks[iflav].mass);

				for(int t=0;t<glbSize[0];t++)
				{
					int ic = iflav*glbSize[0]+t;
					master_fprintf(file2,"%d %+16.16lg %+16.16lg\n",t,corr[ic][RE]*norm,corr[ic][IM]*norm);
				}
				master_fprintf(file2,"\n");
			} 

		}
	close_file(file2);

	nissa_free(corr);

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
