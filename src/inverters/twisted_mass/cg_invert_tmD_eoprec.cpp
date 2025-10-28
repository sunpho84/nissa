#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <optional>
#include <string.h>

#ifdef USE_TMLQCD
# include "base/tmLQCD_bridge.hpp"
#endif

# include "base/multiGridParams.hpp"

#ifdef USE_QUDA
# include "base/quda_bridge.hpp"
#endif

#ifdef USE_DDALPHAAMG
# include "base/DDalphaAMG_bridge.hpp"
#endif

#include "dirac_operators/tmQ/dirac_operator_tmQ.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_mix.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"


#include <dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec_portable.hpp>
#include "cg_64_invert_tmD_eoprec.hpp"

namespace nissa
{
  //Refers to the doc: "doc/eo_inverter.lyx" for explenations
  
  //Prepare the source according to Equation (8.b)
  void inv_tmD_cg_eoprec_prepare_source(OddField<spincolor>& varphi,
					const EoField<quad_su3>& conf_eos,
					const EvnField<spincolor>& eq8a,
					const OddField<spincolor>& source_odd)
  {
    tmn2Deo_or_tmn2Doe_eos(varphi,conf_eos,eq8a);
    
    PAR(0,locVolh,
	CAPTURE(TO_WRITE(varphi),
		TO_READ(source_odd)),
	iEo,
	{
	  for(int id=0;id<NDIRAC/2;id++)
	    for(int ic=0;ic<NCOL;ic++)
	      for(int ri=0;ri<2;ri++)
		{ //gamma5 is explicitly wrote
		  varphi[iEo][id  ][ic][ri]=+source_odd[iEo][id  ][ic][ri]+varphi[iEo][id  ][ic][ri]*0.5;
		  varphi[iEo][id+2][ic][ri]=-source_odd[iEo][id+2][ic][ri]-varphi[iEo][id+2][ic][ri]*0.5;
		}
	});
  }
  
  //invert Koo defined in equation (7)
  void inv_tmDkern_eoprec_square_eos_cg(OddField<spincolor>& sol,
					std::optional<OddField<spincolor>> guess,
					const EoField<quad_su3>& conf,
					const double& kappa,
					const double& mass,
					const int& nitermax,
					const double& residue,
					const OddField<spincolor>& source)
  {
    if(use_128_bit_precision)
      {
	CRASH("reimplement");
	//inv_tmDkern_eoprec_square_eos_cg_128(sol,guess,conf,kappa,mass,nitermax,residue,source);
      }
    else
      inv_tmDkern_eoprec_square_eos_cg_64(sol,guess,conf,kappa,mass,nitermax,residue,source);
  }
  
  //Eq.(11) up to last piece
  void inv_tmD_cg_eoprec_almost_reco_sol(EvnField<spincolor>& varphi,
					 const EoField<quad_su3>& conf_eos,
					 const OddField<spincolor>& sol_odd,
					 const EvnField<spincolor>& source_evn)
  {
    tmn2Deo_or_tmn2Doe_eos(varphi,conf_eos,sol_odd);
    
    PAR(0,locVolh,
	CAPTURE(TO_WRITE(varphi),
		TO_READ(source_evn)),
	iEo,
	{
	  for(int id=0;id<NDIRAC;id++)
	    for(int ic=0;ic<NCOL;ic++)
	      for(int ri=0;ri<2;ri++)
		varphi[iEo][id][ic][ri]=source_evn[iEo][id][ic][ri]+varphi[iEo][id][ic][ri]*0.5;
	});
  }
  
  //Invert twisted mass operator using e/o preconditioning.
  void inv_tmD_cg_eoprec_native(LxField<spincolor>& solution_lx,
				std::optional<OddField<spincolor>> guess_Koo,
				const LxField<quad_su3>& conf_lx,
				const double& kappa,
				const double& mass,
				const int& nitermax,
				const double& residue,
				const LxField<spincolor>& source_lx)
  {
    if(not use_eo_geom) CRASH("eo geometry needed to use cg_eoprec");
    
    //prepare the e/o split version of the source
    EoField<spincolor> source_eos("source_eos",WITH_HALO);
    split_lx_vector_into_eo_parts(source_eos,source_lx);
    
    //prepare the e/o split version of the solution
    EoField<spincolor> solution_eos("solution_eos",WITH_HALO);
    
    //prepare the e/o split version of the conf
    EoField<quad_su3> conf_eos("conf_eos",WITH_HALO);
    split_lx_vector_into_eo_parts(conf_eos,conf_lx);
    
    ///////////////////////////////////// invert with e/o preconditioning ///////////////////////////////////
    
    //Equation (8.a)
    EvnField<spincolor> evnTemp("evnTemp",WITH_HALO);
    inv_tmDee_or_oo_eos(evnTemp,kappa,mass,source_eos.evenPart);
    
    //Prepare the source according to Equation (8.b)
    OddField<spincolor> varphi("varphi",WITH_HALO);
    inv_tmD_cg_eoprec_prepare_source(varphi,conf_eos,evnTemp,source_eos.oddPart);
    
    //Equation (9) using solution_eos[EVN] as temporary vector
    OddField<spincolor> oddTemp("oddTemp",WITH_HALO);
    inv_tmDkern_eoprec_square_eos_cg(oddTemp,guess_Koo,conf_eos,kappa,mass,nitermax,residue,varphi);
    if(guess_Koo) *guess_Koo=oddTemp; //if a guess was passed, return new one
    
    //Equation (10)
    tmDkern_eoprec_eos(solution_eos.oddPart,solution_eos.evenPart,conf_eos,kappa,-mass,oddTemp);
    
    //Equation (11)
    inv_tmD_cg_eoprec_almost_reco_sol(evnTemp,conf_eos,solution_eos.oddPart,source_eos.evenPart);
    inv_tmDee_or_oo_eos(solution_eos.evenPart,kappa,mass,evnTemp);
    
    /////////////////////////// paste the e/o parts of the solution together and free ///////////////////
    
    paste_eo_parts_into_lx_vector(solution_lx,solution_eos);
  }
  
  void inv_tmD_cg_eoprec(LxField<spincolor>& solution_lx,
			 std::optional<OddField<spincolor>> guess_Koo,
			 const LxField<quad_su3>& conf_lx,
			 const double& kappa,
			 const double& mass,
			 const int& nitermax,
			 const double& residue,
			 const LxField<spincolor>& source_lx)
  {
    /// Keep track of convergence
    bool solved=false;
    
    if(multiGrid::checkIfMultiGridAvailableAndRequired(mass) and not solved)
      {
	[[maybe_unused]]
	const double cSW=0;
	
	double call_time=take_time();
#ifdef USE_QUDA
	CRASH("reimplement");
	// solved=quda_iface::solve_tmD(solution_lx,conf_lx,kappa,cSW,mass,nitermax,residue,source_lx);
#elif defined(USE_DDALPHAAMG)
	CRASH("reimplement the handle with DDalpha");
	(void)&cSW;//avoid warning
	//solved=DD::solve(solution_lx,conf_lx,kappa,cSW,mass,residue,source_lx);
#else
	CRASH("How possible!");
#endif
	MASTER_PRINTF("calling multigrid to solve took %lg s\n",take_time()-call_time);
      }

#ifdef USE_TMLQCD
    if(checkIfTmLQCDAvailableAndRequired() and not solved)
      CRASH("Not yet implemented");
#endif
    
    if(not solved)
      inv_tmD_cg_eoprec_native(solution_lx,guess_Koo,conf_lx,kappa,mass,nitermax,residue,source_lx);
    
    //check solution
    double check_time=take_time();
    LxField<spincolor> residueVec("residueVec");
    LxField<spincolor> temp("temp",WITH_HALO);
    temp=solution_lx;
    apply_tmQ(residueVec,conf_lx,kappa,mass,temp);
    PAR(0,locVol,
	CAPTURE(TO_WRITE(residueVec)),
	ivol,
	{
	  safe_dirac_prod_spincolor(residueVec[ivol],base_gamma[5],residueVec[ivol]);
	});
    residueVec-=source_lx;
    
    /// Source L2 norm
    const double sourceNorm2=source_lx.norm2();
    
    /// Residue L2 norm
    const double residueNorm2=residueVec.norm2();
    
    MASTER_PRINTF("check solution, residue: %lg, target one: %lg checked in %lg s\n",residueNorm2/sourceNorm2,residue,take_time()-check_time);
  }
}

