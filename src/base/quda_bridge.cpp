#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_QUDA_BRIDGE
 #include "quda_bridge.hpp"

#include "base/cuda.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "threads/threads.hpp"

namespace quda_iface
{
  using su3_ptr=su3*;
  using quda_conf_t=su3_ptr[NDIM];
  
  using namespace nissa;
  
  /// Lookup tables to map from nissa to quda and vice-versa
  int* loclx_of_quda=nullptr;
  int* quda_of_loclx=nullptr;
  
  /// Conf used to remap
  quda_conf_t quda_conf{};
  
  /// Spincolor used to remap
  spincolor *spincolor_in=nullptr;
  spincolor *spincolor_out=nullptr;
  
  /// Return the rank of the given quda coords
  int get_rank_of_quda_coords(const int *coords,void *fdata)
  {
    /// Coordinates
    int c[NDIM];
    for(int mu=0;mu<NDIM;mu++)
      c[mu]=coords[std::array<int,NDIM>{3,0,1,2}[mu]];
    
    /// Result
    int out;
    MPI_Cart_rank(cart_comm,c,&out);
    
    printf("rank %d, {%d %d %d %d}->%d\n",rank,c[0],c[1],c[2],c[3],out);
    
    return out;
  }
  
  /// Set the verbosity
  void set_quda_verbosity()
  {
    switch(verbosity_lv)
      {
      case 0:
	inv_param.verbosity=QUDA_SILENT;
	break;
      case 1:
	inv_param.verbosity=QUDA_SUMMARIZE;
	break;
      case 2:
	inv_param.verbosity=QUDA_VERBOSE;
	break;
      default:
	inv_param.verbosity=QUDA_DEBUG_VERBOSE;
      }
    
    setVerbosityQuda(QUDA_SUMMARIZE,"# QUDA: ",stdout);
  }
  
  /// Initialize QUDA
  void initialize()
  {
    if(not inited)
      {
	master_printf("Initializing QUDA\n");
	
	if(QUDA_VERSION_MAJOR==0 and QUDA_VERSION_MINOR<7)
	  crash("minimum QUDA version required is 0.7.0");
	
	/////////////////////////////////////////////////////////////////
	
	quda_of_loclx=nissa_malloc("quda_of_loclx",locVol,int);
	loclx_of_quda=nissa_malloc("loclx_of_quda",locVol,int);
	
	for(int ivol=0;ivol<locVol;ivol++)
	  {
	    const coords& c=locCoordOfLoclx[ivol];
	    const coords& l=locSize;
	    
	    int itmp=0;
	    for(int mu=0;mu<NDIM;mu++)
	      {
		const int nu=std::array<int,NDIM>{0,3,2,1}[mu];
		itmp=itmp*l[nu]+c[nu];
	      }
	    const int quda=loclx_parity[ivol]*locVolh+itmp/2;
	    
	    if(quda<0 or quda>=locVol)
	      crash("quda %d remapping to ivol %d not in range [0,%d]",quda,ivol,locVol);
	    
	    quda_of_loclx[ivol]=quda;
	    loclx_of_quda[quda]=ivol;
	  }
	
	for(int mu=0;mu<NDIM;mu++)
	  quda_conf[mu]=nissa_malloc("gauge_cuda",locVol,su3);
	
	spincolor_in=nissa_malloc("spincolor_in",locVol,spincolor);
	spincolor_out=nissa_malloc("spincolor_out",locVol,spincolor);
	
	set_quda_verbosity();
	
	/////////////////////////// gauge params ////////////////////////////////
	
	gauge_param=newQudaGaugeParam();
	
	for(int mu=0;mu<NDIM;mu++)
	  gauge_param.X[mu]=locSize[std::array<int,NDIM>{1,2,3,0}[mu]];
	
	gauge_param.anisotropy=1.0;
	
	gauge_param.type=QUDA_WILSON_LINKS;
	gauge_param.gauge_order=QUDA_QDP_GAUGE_ORDER;
	
	gauge_param.t_boundary=QUDA_PERIODIC_T;
	
	gauge_param.cpu_prec=QUDA_DOUBLE_PRECISION;
	gauge_param.cuda_prec=QUDA_DOUBLE_PRECISION;
	gauge_param.cuda_prec_sloppy=QUDA_SINGLE_PRECISION;
	gauge_param.cuda_prec_precondition=QUDA_HALF_PRECISION;
	gauge_param.cuda_prec_refinement_sloppy=QUDA_SINGLE_PRECISION;
	
	gauge_param.reconstruct=QUDA_RECONSTRUCT_NO;
	gauge_param.reconstruct_sloppy=QUDA_RECONSTRUCT_NO;
	gauge_param.reconstruct_precondition=QUDA_RECONSTRUCT_NO;
	gauge_param.reconstruct_refinement_sloppy=QUDA_RECONSTRUCT_NO;
	
	gauge_param.gauge_fix=QUDA_GAUGE_FIXED_NO;
	
	gauge_param.ga_pad=0;
	for(int mu=0;mu<NDIM;mu++)
	  {
	    int surf_size=locVol/locSize[mu]/2;
	    gauge_param.ga_pad=std::max(gauge_param.ga_pad,surf_size);
	  }
	
	///////////////////////////// inverter parameters ////////////////////////////////////
	
	inv_param=newQudaInvertParam();
	
	inv_param.Ls=1;
	
	inv_param.dagger=QUDA_DAG_NO;
	inv_param.mass_normalization=QUDA_MASS_NORMALIZATION; ///Check
	inv_param.solver_normalization=QUDA_DEFAULT_NORMALIZATION;
	
	inv_param.pipeline=0;
	inv_param.gcrNkrylov=20;
	
	inv_param.residual_type=QUDA_L2_RELATIVE_RESIDUAL;
	inv_param.tol_hq=0.1;
	inv_param.reliable_delta=1e-3; // ignored by multi-shift solver
	inv_param.use_sloppy_partial_accumulator=0;
	
	// domain decomposition preconditioner parameters
	inv_param.inv_type_precondition=QUDA_CG_INVERTER;
	inv_param.schwarz_type=QUDA_ADDITIVE_SCHWARZ;
	inv_param.precondition_cycle=1;
	inv_param.tol_precondition=0.1;
	inv_param.maxiter_precondition=10;
	inv_param.verbosity_precondition=QUDA_SILENT;
	if(verbosity_lv>=2)
	  inv_param.verbosity_precondition=QUDA_VERBOSE;
	
	inv_param.omega=1.0;
	
	inv_param.cpu_prec=QUDA_DOUBLE_PRECISION;
	inv_param.cuda_prec=QUDA_DOUBLE_PRECISION;
	inv_param.cuda_prec_sloppy=QUDA_SINGLE_PRECISION;
	inv_param.cuda_prec_refinement_sloppy=QUDA_SINGLE_PRECISION;
	inv_param.cuda_prec_precondition=QUDA_HALF_PRECISION;
	
	inv_param.clover_cpu_prec=QUDA_DOUBLE_PRECISION;
	inv_param.clover_cuda_prec=QUDA_DOUBLE_PRECISION;
	inv_param.clover_cuda_prec_sloppy=QUDA_SINGLE_PRECISION;
	inv_param.clover_cuda_prec_precondition=QUDA_HALF_PRECISION;
	inv_param.clover_cuda_prec_refinement_sloppy=QUDA_SINGLE_PRECISION;
	
	inv_param.preserve_source=QUDA_PRESERVE_SOURCE_YES;
	inv_param.gamma_basis=QUDA_CHIRAL_GAMMA_BASIS;
	inv_param.dirac_order=QUDA_DIRAC_ORDER;
	
	inv_param.input_location=QUDA_CPU_FIELD_LOCATION;
	inv_param.output_location=QUDA_CPU_FIELD_LOCATION;
	
	inv_param.tune=QUDA_TUNE_YES;
	
	inv_param.sp_pad=0;
	inv_param.cl_pad=0;
	
	inv_mg_param=newQudaInvertParam();
	quda_mg_param=newQudaMultigridParam();
	
	int grid[NDIM];
	for(int mu=0;mu<NDIM;mu++)
	  grid[mu]=nrank_dir[std::array<int,NDIM>{1,2,3,0}[mu]];
	
	initCommsGridQuda(NDIM,grid,get_rank_of_quda_coords,NULL);
	
	initQuda(iCudaDevice);
	
	inited=1;
      }
  }
  
  /// Finalize QUDA
  void finalize()
  {
    if(inited)
      {
	master_printf("Finalizing QUDA\n");
	
	nissa_free(loclx_of_quda);
	nissa_free(quda_of_loclx);
	
	// destroys the preconditioner if it was created
	if(quda_mg_preconditioner!=NULL)
	  {
	    destroyMultigridQuda(quda_mg_preconditioner);
	    quda_mg_preconditioner=NULL;
	  }
	
	for(int mu=0;mu<NDIM;mu++)
	  nissa_free(quda_conf[mu]);
	
	nissa_free(spincolor_in);
	nissa_free(spincolor_out);
	
	//free the clover and gauge conf
	freeGaugeQuda();
	freeCloverQuda();
	endQuda();
	
	inited=false;
      }
  }
  
  /// Loads the clover term
  void load_clover_term(QudaInvertParam* inv_param)
  {
    freeCloverQuda();
    loadCloverQuda(NULL,NULL,inv_param);
  }
  
  /// Reorder conf into QUDA format
  void remap_nissa_to_quda(quda_conf_t out,quad_su3 *in)
  {
    master_printf("%s\n",__FUNCTION__);
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	const int iquda=quda_of_loclx[ivol];
	
	for(int nu=0;nu<NDIM;nu++)
	  su3_copy(out[std::array<int,NDIM>{3,0,1,2}[nu]][iquda],in[ivol][nu]);
      }
    NISSA_PARALLEL_LOOP_END;
    
    THREAD_BARRIER();
  }
  
  /// Reorder spinor to QUDA format
  void remap_nissa_to_quda(spincolor *out,spincolor *in)
  {
    master_printf("%s\n",__FUNCTION__);
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	const int iquda=quda_of_loclx[ivol];
	spincolor_copy(out[iquda],in[ivol]);
      }
    NISSA_PARALLEL_LOOP_END;
    
    THREAD_BARRIER();
  }
  
  /// Reorder spinor from QUDA format
  void remap_quda_to_nissa(spincolor *out,spincolor *in)
  {
    NISSA_PARALLEL_LOOP(iquda,0,locVol)
      {
	const int ivol=loclx_of_quda[iquda];
	spincolor_copy(out[ivol],in[iquda]);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  /// Load a gauge conf
  void load_conf(quad_su3 *nissa_conf)
  {
    master_printf("freeing the QUDA gauge conf\n");
    freeGaugeQuda();
    
    remap_nissa_to_quda(quda_conf,nissa_conf);
    master_printf("loading to QUDA the gauge conf\n");
    loadGaugeQuda((void*)quda_conf,&gauge_param);
    
    double plaq;
    plaqQuda(&plaq);
    master_printf("loaded, plaquette: %.16lg\n",plaq);
  }
  
  /// Sets the sloppy precision
  void set_sloppy_prec(const QudaPrecision sloppy_prec)
  {
    master_printf("Quda sloppy precision: ");
    
    switch(sloppy_prec)
      {
      case QUDA_DOUBLE_PRECISION:
	master_printf("double");
	break;
      case QUDA_SINGLE_PRECISION:
	master_printf("single");
	break;
      case QUDA_HALF_PRECISION:
	master_printf("half");
	break;
      case QUDA_QUARTER_PRECISION:
	master_printf("quarter");
	break;
      case QUDA_INVALID_PRECISION:
	crash("invalid precision");
      };
    master_printf("\n");
  
  gauge_param.cuda_prec_sloppy=
    gauge_param.cuda_prec_refinement_sloppy=
    inv_param.cuda_prec_sloppy=
    inv_param.clover_cuda_prec_sloppy=
    inv_param.clover_cuda_prec_refinement_sloppy=
    sloppy_prec;
  }
  
  /// Apply the dirac operator
  void apply_tmD(spincolor *out,quad_su3 *conf,double kappa,double mu,spincolor *in)
  {
    load_conf(conf);
    
    master_printf("setting pars\n");
    
    inv_param.kappa=kappa;
    inv_param.dslash_type=QUDA_TWISTED_MASS_DSLASH;
    inv_param.solution_type=QUDA_MAT_SOLUTION;
    
    //minus due to different gamma5 definition
    inv_param.mu=-mu;
    inv_param.epsilon=0.0;
    
    inv_param.twist_flavor=QUDA_TWIST_SINGLET;
    inv_param.Ls=1;
    
    remap_nissa_to_quda(spincolor_in,in);
    MatQuda(spincolor_out,spincolor_in,&inv_param);
    remap_quda_to_nissa(out,spincolor_out);
  }
  
  void solve(spincolor *sol,quad_su3 *conf,double kappa,double mu,int niter,double residue,spincolor *source)
  {
    load_conf(conf);
    
    inv_param.kappa=kappa;
    
    inv_param.dslash_type=QUDA_TWISTED_MASS_DSLASH;
    inv_param.matpc_type=QUDA_MATPC_EVEN_EVEN_ASYMMETRIC;
    inv_param.solution_type=QUDA_MAT_SOLUTION;
    
    inv_param.inv_type=QUDA_CG_INVERTER;
    inv_param.solve_type=QUDA_NORMERR_PC_SOLVE;
    
    //minus due to different gamma5 definition
    inv_param.mu=-mu;// /(2.0*kappa); /// Check kappa
    inv_param.epsilon=0.0;
    
    inv_param.twist_flavor=QUDA_TWIST_SINGLET;
    master_printf("Residue: %lg\n",residue);
    inv_param.tol=sqrt(residue);
    inv_param.maxiter=niter;
    inv_param.Ls=1;
    
    set_quda_verbosity();
    
    remap_nissa_to_quda(spincolor_in,source);
    
    invertQuda(spincolor_out,spincolor_in,&inv_param);
    
    master_printf("# QUDA solved in: %i iter / %g secs = %g Gflops\n",inv_param.iter,inv_param.secs,inv_param.gflops/inv_param.secs);
    
    remap_quda_to_nissa(sol,spincolor_out);
  }
}
