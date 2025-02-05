#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#ifdef USE_TMLQCD
# include "base/tmLQCD_bridge.hpp"
#endif

#ifdef USE_QUDA
# include "base/quda_bridge.hpp"
#endif

#ifdef USE_DDALPHAAMG
# include "base/DDalphaAMG_bridge.hpp"
#endif

#include "dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec.hpp"
#include "dirac_operators/tmclovQ/dirac_operator_tmclovQ.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_mix.hpp"
#include "inverters/twisted_mass/cg_invert_tmD_eoprec.hpp"
#include "new_types/su3_op.hpp"
#include "operations/su3_paths/clover_term.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

#include "cg_64_invert_tmclovD_eoprec.hpp"

namespace nissa
{
  //Refers to the tmD_eoprec
  
  //invert Koo defined in equation (7)
  void inv_tmclovDkern_eoprec_square_eos_cg(OddField<spincolor>& sol,
					    std::optional<OddField<spincolor>> guess,
					    const EoField<quad_su3>& conf,
					    const double& kappa,
					    const double& cSW,
					    const OddField<clover_term_t>& Cl_odd,
					    const EvnField<inv_clover_term_t>& invCl_evn,
					    const double& mass,
					    const int& nitermax,
					    const double& residue,
					    const OddField<spincolor>& source)
  {
    if(use_128_bit_precision)
      {
	CRASH("reimplement");
	//inv_tmclovDkern_eoprec_square_eos_cg_128(sol,guess,conf,kappa,cSW,Cl_odd,invCl_evn,mass,nitermax,residue,source);
      }
    else inv_tmclovDkern_eoprec_square_eos_cg_64(sol,guess,conf,kappa,Cl_odd,invCl_evn,mass,nitermax,residue,source);
  }
  
  //Invert twisted clover operator using e/o preconditioning.
  void inv_tmclovD_cg_eoprec_native(LxField<spincolor>& solution_lx,
				    std::optional<OddField<spincolor>> guess_Koo,
				    const LxField<quad_su3>& conf_lx,
				    const double& kappa,
				    const double& cSW,
				    const LxField<clover_term_t>& Cl_lx,
				    const LxField<inv_clover_term_t>* ext_invCl_lx,
				    const double& mass,
				    const int& nitermax,
				    const double& residue,
				    const LxField<spincolor>& source_lx)
  {
    // set_borders_invalid(conf_lx);
    // communicate_lx_quad_su3_borders(conf_lx);
    if(not use_eo_geom)
      CRASH("eo geometry needed to use cg_eoprec");
    
    const LxField<inv_clover_term_t> *_invCl_lx;
    if(ext_invCl_lx)
      _invCl_lx=ext_invCl_lx;
    else
      {
	LxField<inv_clover_term_t>* tmp=new LxField<inv_clover_term_t>("invCl");
	_invCl_lx=tmp;
	invert_twisted_clover_term(*tmp,mass,kappa,Cl_lx);
      }
    const LxField<inv_clover_term_t> invCl_lx=*_invCl_lx;
    
    //prepare the e/o split version of the source
    EoField<spincolor> source_eos("source_eos",WITH_HALO);
    split_lx_vector_into_eo_parts(source_eos,source_lx);
    
    //prepare the e/o split version of the solution
    EoField<spincolor> solution_eos("solution_eos",WITH_HALO);
    
    //prepare the e/o split version of the conf
    EoField<quad_su3> conf_eos("conf_eos",WITH_HALO);
    split_lx_vector_into_eo_parts(conf_eos,conf_lx);
    
    //prepare the e/o split version of the clover term
    OddField<clover_term_t> Cl_odd("Cl_odd");
    get_evn_or_odd_part_of_lx_vector(Cl_odd,Cl_lx,ODD_SITES);
    
    //prepare the e/o split version of the clover term
    EvnField<inv_clover_term_t> invCl_evn("invCl_evn");
    get_evn_or_odd_part_of_lx_vector(invCl_evn,invCl_lx,EVEN_SITES);
    
    ///////////////////////////////////// invert with e/o preconditioning ///////////////////////////////////
    
    //Equation (8.a)
    EvnField<spincolor> evnTemp("evnTemp",WITH_HALO);
    inv_tmclovDee_or_oo_eos(evnTemp,invCl_evn,false,source_eos.evenPart);
    
    //Equation (8.b)
    OddField<spincolor> varphi("varphi",WITH_HALO);
    inv_tmD_cg_eoprec_prepare_source(varphi,conf_eos,evnTemp,source_eos.oddPart);
    
    //Equation (9) using solution_eos[EVN] as temporary vector
    OddField<spincolor> oddTemp("oddTemp",WITH_HALO);
    inv_tmclovDkern_eoprec_square_eos_cg(oddTemp,guess_Koo,conf_eos,kappa,cSW,Cl_odd,invCl_evn,mass,nitermax,residue,varphi);
    if(guess_Koo) (*guess_Koo)=oddTemp; //if a guess was passed, return new one
    
    //Equation (10)
    tmclovDkern_eoprec_eos(solution_eos.oddPart,solution_eos.evenPart,conf_eos,kappa,Cl_odd,invCl_evn,true,mass,oddTemp);
    
    //Equation (11)
    inv_tmD_cg_eoprec_almost_reco_sol(evnTemp,conf_eos,solution_eos.oddPart,source_eos.evenPart);
    inv_tmclovDee_or_oo_eos(solution_eos.evenPart,invCl_evn,false,evnTemp);
    
    /////////////////////////// paste the e/o parts of the solution together and free ///////////////////
    
    paste_eo_parts_into_lx_vector(solution_lx,solution_eos);
    
    // //check for residual
    // LxField<spincolor> check_res("check_res",WITH_HALO);
    
    // //multiply by g5*D
    // apply_tmclovQ(check_res,conf_lx,kappa,Cl_lx,mass,solution_lx);
    
    // //remove g5 and take the difference with source
    // PAR(0,locVol,
    // 	CAPTURE(TO_WRITE(check_res),
    // 		TO_READ(source_lx)),ivol,
    //   {
    // 	const double mg5[2]={-1,1};
    // 	for(int high_low=0;high_low<2;high_low++)
    // 	  for(int id=high_low*NDIRAC/2;id<(high_low+1)*NDIRAC/2;id++)
    // 	    color_summ_the_prod_double(check_res[ivol][id],source_lx[ivol][id],mg5[high_low]);
    //   });
    
    // //compute residual and print
    // const double real_residue=check_res.norm2();
    // if(real_residue>residue*1.1)
    //   WARNING("preconditioned tmclovD solver, asked for residual: %lg, obtained %lg\n\n",residue,real_residue);
    
    if(ext_invCl_lx==nullptr)
      delete _invCl_lx;
  }
  
  void inv_tmclovD_cg_eoprec(LxField<spincolor>& solution_lx,
			     std::optional<OddField<spincolor>> guess_Koo,
			     const LxField<quad_su3>& conf_lx,
			     const double& kappa,
			     const LxField<clover_term_t>& Cl_lx,
			     const LxField<inv_clover_term_t>* ext_invCl_lx,
			     const double& cSW,
			     const double& mass,
			     const int& nitermax,
			     const double& targResidue,
			     const LxField<spincolor>& source_lx)
  {
    /// Keep track of convergence
    bool solved=false;
    
#ifdef USE_QUDA
    if(use_quda and not solved)
      {
	const double call_time=take_time();
	solved=quda_iface::solve_tmD(solution_lx,conf_lx,kappa,cSW,mass,nitermax,targResidue,source_lx);
	MASTER_PRINTF("calling quda to solve took %lg s\n",take_time()-call_time);
      }
#endif
    
#ifdef USE_DDALPHAAMG
    if(multiGrid::checkIfMultiGridAvailableAndRequired(mass) and not solved)
      {
	if(source_lx.spaceTimeLayout!=SpaceTimeLayout::CPU)
	  CRASH("wrong layout");
	
	const double call_time=take_time();
	solved=DD::solve((spincolor*)solution_lx.getPtr<MemorySpace::CPU>(),
			 (quad_su3*)conf_lx.getPtr<MemorySpace::CPU>(),
			 kappa,cSW,mass,targResidue,(spincolor*)source_lx.getPtr<MemorySpace::CPU>());
	MASTER_PRINTF("calling DDalphaAMG to solve took %lg s\n",take_time()-call_time);
      }
#endif
    
    if(not solved)
      inv_tmclovD_cg_eoprec_native(solution_lx,guess_Koo,conf_lx,kappa,cSW,Cl_lx,ext_invCl_lx,mass,nitermax,targResidue,source_lx);
    
    if(check_inversion_residue)
      {
	double check_time=take_time();
	LxField<spincolor> residueVec("residueVec");
	apply_tmclovQ(residueVec,conf_lx,kappa,Cl_lx,mass,solution_lx);
	
	PAR(0,locVol,
	    CAPTURE(TO_WRITE(residueVec)),ivol,
	    {
	      unsafe_dirac_prod_spincolor(residueVec[ivol],base_gamma[5],residueVec[ivol]);
	    });
	
	residueVec-=source_lx;
	
	/// Residue L2 norm
	const double residueNorm2=
	  residueVec.norm2();
	
	/// Source L2 norm
	const double sourceNorm2=
	  source_lx.norm2();
	
	const double residue=
	  residueNorm2/sourceNorm2;
	
	MASTER_PRINTF(" Check solution, source norm2: %lg, residue: %lg, target one: %lg checked in %lg s\n",
		      sourceNorm2,residue,targResidue,take_time()-check_time);
	
	double ref;
	if(targResidue<1e-29)
	  {
	    ref=pow(10,-inversion_residue_heavy_qualify_odg)*sqrt(glbVol);
	    MASTER_PRINTF("  Heavy-quark residue, comparing with absolute threshold %lg\n",ref);
	  }
	else
	  {
	    ref=targResidue*pow(10,inversion_residue_threshold_odg);
	    MASTER_PRINTF("  Comparing with %d orders larger residue: %lg\n",inversion_residue_threshold_odg,ref);
	  }
	
	if(residue>ref)
	  {
	    const char* txt=
	      "residue is larger than threshold";
	    
	    if(check_inversion_residue>1)
	      CRASH("%s",txt);
	    else
	      WARNING("%s",txt);
	  }
	else
 	  MASTER_PRINTF(GREEN_HIGHLIGHT "    Inversion passed the residue check " DO_NOT_HIGHLIGHT "\n");
      }
  }
}
