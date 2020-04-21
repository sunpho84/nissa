#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_QUDA_BRIDGE
 #include "quda_bridge.hpp"

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"

namespace quda_iface
{
  /// Direction to be used to refer to the nissa one
  CUDA_MANAGED int nissa_dir_of_quda[NDIM]={1,2,3,0};
  CUDA_MANAGED int quda_dir_of_nissa[NDIM]={3,0,1,2};
  
  /// Return the rank of the given quda coords
  int get_rank_of_quda_coords(const int *coords,void *fdata)
  {
    /// Coordinates
    int c[NDIM];
    for(int mu=0;mu<NDIM;mu++)
      c[mu]=nissa::rank_coord[nissa_dir_of_quda[mu]];
    
    /// Result
    int out;
    MPI_Cart_rank(nissa::cart_comm,c,&out);
    
    return out;
  }
  
  void initialize()
  {
    if(not inited)
      {
	if(QUDA_VERSION_MAJOR==0 and QUDA_VERSION_MINOR<7)
	  crash("minimum QUDA version required is 0.7.0");
	
	////////////////////////////// verbosity ///////////////////////////////////
	
	switch(nissa::verbosity_lv)
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
	
	/////////////////////////// gauge params ////////////////////////////////
	
	gauge_param=newQudaGaugeParam();
	
	for(int mu=0;mu<NDIM;mu++)
	  gauge_param.X[mu]=nissa::loc_size[nissa_dir_of_quda[mu]];
	
	gauge_param.anisotropy=1.0;
	gauge_param.type=QUDA_WILSON_LINKS;
	gauge_param.gauge_order=QUDA_QDP_GAUGE_ORDER;
	
	gauge_param.cpu_prec=QUDA_DOUBLE_PRECISION;
	gauge_param.cuda_prec=QUDA_DOUBLE_PRECISION;
	gauge_param.reconstruct=QUDA_RECONSTRUCT_NO;
	gauge_param.cuda_prec_sloppy=QUDA_SINGLE_PRECISION;
	gauge_param.reconstruct_sloppy=QUDA_RECONSTRUCT_NO;
	gauge_param.cuda_prec_precondition=QUDA_HALF_PRECISION;
	gauge_param.reconstruct_precondition=QUDA_RECONSTRUCT_NO;
	gauge_param.cuda_prec_refinement_sloppy=QUDA_SINGLE_PRECISION;
	gauge_param.reconstruct_refinement_sloppy=QUDA_RECONSTRUCT_NO;
	gauge_param.gauge_fix=QUDA_GAUGE_FIXED_NO;
	
	gauge_param.ga_pad=0;
	for(int mu=0;mu<NDIM;mu++)
	  gauge_param.ga_pad=std::max(gauge_param.ga_pad,nissa::bord_dir_vol[mu]/2);
	
	///////////////////////////// inverter parameters ////////////////////////////////////
	
	inv_param=newQudaInvertParam();
	
	inv_param.Ls=1;
	
	inv_param.dagger=QUDA_DAG_NO;
	inv_param.mass_normalization=QUDA_KAPPA_NORMALIZATION;
	inv_param.solver_normalization=QUDA_DEFAULT_NORMALIZATION;
	
	inv_param.pipeline=0;
	inv_param.gcrNkrylov=20;
	
	inv_param.residual_type=(QudaResidualType)QUDA_L2_RELATIVE_RESIDUAL;
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
	if(nissa::verbosity_lv>=2)
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
	  grid[mu]=nissa::nrank_dir[nissa_dir_of_quda[mu]];
	
	initCommsGridQuda(NDIM,grid,get_rank_of_quda_coords,NULL);
	
	int idevice=0;
	master_printf("Fix this\n");
	initQuda(idevice);
	
	inited=1;
      }
  }
  
  void finalize()
  {
    if(inited)
      {
	// destroys the preconditioner if it was created
	if(quda_mg_preconditioner!=NULL)
	  {
	    destroyMultigridQuda(quda_mg_preconditioner);
	    quda_mg_preconditioner=NULL;
	  }
	
	//free the clover and gauge conf
	freeGaugeQuda();
	freeCloverQuda();
	endQuda();
      }
  }
  
  void load_clover_term(QudaInvertParam* inv_param)
  {
    freeCloverQuda();
    loadCloverQuda(NULL,NULL,inv_param);
  }
  
  void load_conf(nissa::quad_su3 *nissa_conf)
  {
    using namespace nissa;
    
    freeGaugeQuda();
    
    //Allocate the temporary buffer
    double *quda_conf=nissa_malloc("gauge_cuda",loc_vol*sizeof(quad_su3)/sizeof(double),double);
    
    GET_THREAD_ID();
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	const coords& c=loc_coord_of_loclx[ivol];
	const coords& l=loc_size;
	const int itmp=c[1]+l[1]*(c[2]+l[2]*(c[3]+l[3]*c[0]));
	const int par=(c[0]+c[1]+c[2]+c[3])&1;
	const int iquda=par*loc_volh+itmp/2;
	
	for(int mu=0;mu<NDIM;mu++)
	  memcpy(quda_conf+(iquda+loc_vol*quda_dir_of_nissa[mu])*NCOL*NCOL*2,nissa_conf[ivol][mu],sizeof(su3));
      }
    NISSA_PARALLEL_LOOP_END;
    
    THREAD_BARRIER();
    
    loadGaugeQuda((void*)quda_conf,&gauge_param);
    
    nissa_free(quda_conf);
  }
}
