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
	crash("reimplement");
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
      crash("eo geometry needed to use cg_eoprec");
    
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
    
    //check for residual
    LxField<spincolor> check_res("check_res",WITH_HALO);
    
    //multiply by g5*D
    apply_tmclovQ(check_res,conf_lx,kappa,Cl_lx,mass,solution_lx);
    
    //remove g5 and take the difference with source
    PAR(0,locVol,
	CAPTURE(TO_WRITE(check_res),
		TO_READ(source_lx)),ivol,
      {
	const double mg5[2]={-1,1};
	for(int high_low=0;high_low<2;high_low++)
	  for(int id=high_low*NDIRAC/2;id<(high_low+1)*NDIRAC/2;id++)
	    color_summ_the_prod_double(check_res[ivol][id],source_lx[ivol][id],mg5[high_low]);
      });
    
    //compute residual and print
    const double real_residue=check_res.norm2();
    if(real_residue>residue*1.1)
      master_printf("WARNING preconditioned tmclovD solver, asked for residual: %lg, obtained %lg\n\n",residue,real_residue);
    
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
			     const double& residue,
			     const LxField<spincolor>& source_lx)
  {
    //memset(send_buf,0,send_buf_size);
    //memset(recv_buf,0,recv_buf_size);
    
  //   {
  //   int ivolIncr,rankIncr;
  //   get_loclx_and_rank_of_coord(ivolIncr,rankIncr,{glbSize[0]-1,8,23,7});
  //   if(rank==rankIncr)
  //     {
  // 	printf("at the beginning of the solver, the bulk of the conf on rank %d is:\n",rank);
  // 	su3_print(conf_lx[ivolIncr][0]);
  //     }
    
  // }
    
  //   const int incrSite=533223;
    
  //   master_printf("conf_lx pointer: %p\n",conf_lx);
    
  //   if(rank==0)
  //     {
  // 	printf("at the beginning there was\n");
  // 	su3_print(conf_lx[incrSite][0]);
  // 	su3_put_to_id(conf_lx[incrSite][0]);
  // 	printf("darkness\n");
  // 	su3_print(conf_lx[incrSite][0]);
  //     }
  //   set_borders_invalid(conf_lx);
  //   communicate_lx_quad_su3_borders(conf_lx);
    
  //   if(rank==0)
  //     {
  // 	printf("then, light came\n");
  // 	su3_print(conf_lx[incrSite][0]);
  //     }
    
  //   {
  //     checksum check;
  //     checksum_compute_nissa_data(check,Cl_lx,64,sizeof(clover_term_t));
  //     master_printf("initial checksum of the clover %x %x\n",check[0],check[1]);
  //   }
    
    /// Keep track of convergence
    bool solved=false;
    
#ifdef USE_QUDA
    if(use_quda and not solved)
      {
	const double call_time=take_time();
	solved=quda_iface::solve_tmD(solution_lx,conf_lx,kappa,cSW,mass,nitermax,residue,source_lx);
	master_printf("calling quda to solve took %lg s\n",take_time()-call_time);
      }
#endif
    
#ifdef USE_DDALPHAAMG
    if(multiGrid::checkIfMultiGridAvailableAndRequired(mass) and not solved)
      {
	if(source_lx.fieldLayout!=FieldLayout::CPU)
	  crash("wrong layout");
	
	const double call_time=take_time();
	solved=DD::solve((spincolor*)solution_lx._data,(quad_su3*)conf_lx._data,kappa,cSW,mass,residue,(spincolor*)source_lx._data);
	master_printf("calling DDalphaAMG to solve took %lg s\n",take_time()-call_time);
      }
#endif
    
    // if(checkIfTmLQCDAvailableAndRequired() and not solved)
    //   {
    // 	// //open the file
    // 	// std::string temp_path="temp_input";
    // 	// FILE *ftemp;
    // 	// master_get_temp_file(ftemp,temp_path);
	
    // 	// //preprare the parameters
    // 	// int argc=3,verbose=1,external_id=0;
    // 	// char *argv[3]={(char*)malloc(10),(char*)malloc(10),(char*)malloc(temp_path.length()+1)};
    // 	// sprintf(argv[0],"-");
    // 	// sprintf(argv[1],"-f");
    // 	// sprintf(argv[2],"%s",temp_path.c_str());
	
    // 	//prepare the input file
	
    // 	//initializing
    // 	FILE *ftemp=open_prepare_input_file_for_tmLQCD();
    // 	master_fprintf(ftemp,"\n");
    // 	master_fprintf(ftemp,"2kappamu = %lg\n",2*kappa*mass);
    // 	master_fprintf(ftemp,"kappa = %lg\n",kappa);
    // 	master_fprintf(ftemp,"DebugLevel = 3\n");
    // 	master_fprintf(ftemp,"csw = %lg\n",cSW);
    // 	master_fprintf(ftemp,"\n");
    // 	//master_fprintf(ftemp,"BeginOperator TMWILSON\n");
    // 	master_fprintf(ftemp,"BeginOperator CLOVER\n");
    // 	master_fprintf(ftemp,"  2kappamu = %lg\n",2*kappa*mass);
    // 	master_fprintf(ftemp,"  kappa = %lg\n",kappa);
    // 	master_fprintf(ftemp,"  cSW = %lg\n",cSW);
    // 	master_fprintf(ftemp,"  UseEvenOdd = yes\n");
    // 	master_fprintf(ftemp,"  Solver = CG\n");
    // 	master_fprintf(ftemp,"  SolverPrecision = %lg\n",residue);
    // 	master_fprintf(ftemp,"  MaxSolverIterations = %d\n",nitermax);
    // 	master_fprintf(ftemp,"  AddDownPropagator = no\n");
    // 	master_fprintf(ftemp,"EndOperator\n");
	
    // 	close_file(ftemp);
	
    // 	tmLQCD_init();
    // 	export_gauge_conf_to_tmLQCD(conf_lx);
    // 	tmLQCD::tmLQCD_invert((double*)solution_lx,(double*)(source_lx),0,0);
	
    // 	tmLQCD_finalise();
	
    // 	// for(int i=0;i<argc;i++) free(argv[i]);
	
    // 	//if(rank==0) unlink("invert.input");
	
    // 	// //close the temporary file and remove it
    // 	// if(rank==0)
    // 	//   {
    // 	// 	fclose(ftemp);
    // 	// 	unlink(temp_path.c_str());
    // 	//   }
    // 	crash("Not yet implemented");
    //   }
    
    if(not solved)
      inv_tmclovD_cg_eoprec_native(solution_lx,guess_Koo,conf_lx,kappa,cSW,Cl_lx,ext_invCl_lx,mass,nitermax,residue,source_lx);
    
    if(check_inversion_residue)
      {
    // THREAD_BARRIER();
	//check solution
	double check_time=take_time();
	LxField<spincolor> residueVec("residueVec");
	// checksum check;
	// checksum_compute_nissa_data(check,solution_lx,64,sizeof(spincolor));
	// master_printf("checksum of the solution %x %x\n",check[0],check[1]);
	// checksum_compute_nissa_data(check,conf_lx,64,sizeof(quad_su3));
	// master_printf("checksum of the conf %x %x\n",check[0],check[1]);
	
	// const double sou=source_lx[0][0][0][0];
	// const double sol=solution_lx[0][0][0][0];
	// set_borders_invalid(solution_lx);
	apply_tmclovQ(residueVec,conf_lx,kappa,Cl_lx,mass,solution_lx);
	// const double sola=solution_lx[0][0][0][0];
	// const double soll=solution_lx[locVol][0][0][0];
	// checksum_compute_nissa_data(check,Cl_lx,64,sizeof(clover_term_t));
	// master_printf("checksum of the clover %x %x\n",check[0],check[1]);
	// checksum_compute_nissa_data(check,residueVec,64,sizeof(spincolor));
	// master_printf("checksum of the residue %x %x\n",check[0],check[1]);
	// checksum_compute_nissa_data(check,solution_lx+bord_vol,64,sizeof(spincolor));
	// master_printf("checksum of the solution shifted by bord %x %x\n",check[0],check[1]);
	// const double res=residueVec[0][0][0][0];
	// const double res1=residueVec[loclx_of_coord_list(0,8,23,7)][0][0][0];
	PAR(0,locVol,
	    CAPTURE(TO_WRITE(residueVec)),ivol,
	    {
	      unsafe_dirac_prod_spincolor(residueVec[ivol],base_gamma[5],residueVec[ivol]);
	    });
	
	// const double res5=residueVec[0][0][0][0];
	residueVec-=source_lx;
	// const double ress=residueVec[0][0][0][0];
	// const double resn=spincolor_norm2(residueVec[0]);
	// double resnt=0;
	// bool found=false;
	// for(int ivol=0;ivol<locVol;ivol++)
	//   {
	//     const double contr=spincolor_norm2(residueVec[ivol]);
	//     if(found==false and contr>1e-4)
	//       {
	// 	found=true;
	// 	const coords_t c=locCoordOfLoclx[ivol];
	// 	// master_
	// 	printf("rank %d found first exceeding, %lg at %d %d %d %d\n",rank,contr,c[0],c[1],c[2],c[3]);
	//       }
	//     // resnt+=contr;
	//   }
	// // double resntg;
	// // MPI_Allreduce(&resnt,&resntg,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	
	/// Residue L2 norm
	const double residueNorm2=residueVec.norm2();
	
	/// Source L2 norm
	const double sourceNorm2=source_lx.norm2();
	
	master_printf("check solution, source norm2: %lg, residue: %lg, target one: %lg checked in %lg s\n",
		      sourceNorm2,residueNorm2/sourceNorm2,residue,take_time()-check_time);
	// printf("check rank %d %lg %lg %lg %lg %lg %lg %lg %lg     %lg %lg %lg\n",rank,sou,sol,sola,soll,res,res1,res5,ress,resn,resnt,resntg);
	
    // if(rank==0)
    //   {
    // 	printf("now\n");
    // 	su3_print(conf_lx[incrSite][0]);
    //   }
	
      }
  }
}
