#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_QUDA_BRIDGE
# include "quda_bridge.hpp"

#include "base/cuda.hpp"
#include "base/export_conf_to_external_lib.hpp"
#include "base/multiGridParams.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "threads/threads.hpp"

#                      include "hmc/backfield.hpp"

#ifdef ENABLE_QUDA_BYPASS

extern "C"
{
#define QUDA_BYPASS(RET,ARGS...)			\
  ARGS							\
  {							\
    master_printf("Faking %s\n",__PRETTY_FUNCTION__);\
    						     \
    return RET;					     \
  }
  
  QUDA_BYPASS(,void setVerbosityQuda(QudaVerbosity verbosity,const char prefix[],FILE *outfile))
  QUDA_BYPASS({},QudaGaugeParam newQudaGaugeParam(void))
  QUDA_BYPASS({},QudaInvertParam newQudaInvertParam(void))
  QUDA_BYPASS({},QudaMultigridParam newQudaMultigridParam(void))
  QUDA_BYPASS({},QudaEigParam newQudaEigParam(void))
  QUDA_BYPASS(,void initCommsGridQuda(int nDim, const int *dims, QudaCommsMap func, void *fdata))
  QUDA_BYPASS(,void initQuda(int device))
  QUDA_BYPASS(,void destroyMultigridQuda(void *mg_instance))
  QUDA_BYPASS(,void freeGaugeQuda(void))
  QUDA_BYPASS(,void freeCloverQuda(void))
  QUDA_BYPASS(,void endQuda(void))
  QUDA_BYPASS(,void MatQuda(void *h_out, void *h_in, QudaInvertParam *inv_param))
  QUDA_BYPASS({},void* newMultigridQuda(QudaMultigridParam *param))
  QUDA_BYPASS(,void loadCloverQuda(void *h_clover, void *h_clovinv,QudaInvertParam *inv_param))
  QUDA_BYPASS(,void invertQuda(void *h_x, void *h_b, QudaInvertParam *param))
  QUDA_BYPASS(,void printQudaInvertParam(QudaInvertParam *param))
  QUDA_BYPASS(,void loadGaugeQuda(void *h_gauge, QudaGaugeParam *param))
  QUDA_BYPASS(,void plaqQuda(double plaq[3]))
}

namespace nissa
{
  int iCudaDevice;
}

#endif

#ifndef USE_CUDA
namespace nissa
{
  int iCudaDevice;
}
#endif

namespace quda_iface
{
  using namespace nissa;
  
  /// Lookup tables to map from nissa to quda and vice-versa
  CUDA_MANAGED int* loclx_of_quda=nullptr;
  CUDA_MANAGED int* quda_of_loclx=nullptr;
  
  /// Color used to remap
  color *color_in=nullptr;
  color *color_out=nullptr;
  
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
  QudaVerbosity get_quda_verbosity()
  {
    switch(verbosity_lv)
      {
      case 0:
	return QUDA_SILENT;
	break;
      case 1:
	return QUDA_SUMMARIZE;
	break;
      case 2:
	return QUDA_VERBOSE;
	break;
      default:
	return QUDA_DEBUG_VERBOSE;
      }
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
	    const coords_t& c=locCoordOfLoclx[ivol];
	    const coords_t& l=locSize;
	    
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
	  quda_conf[mu]=nissa_malloc("quda_conf",locVol,su3);
	master_printf("allocating quda_conf\n");
	
	spincolor_in=nissa_malloc("spincolor_in",locVol,spincolor);
	spincolor_out=nissa_malloc("spincolor_out",locVol,spincolor);
	
	color_in=nissa_malloc("color_in",locVol,color);
	color_out=nissa_malloc("color_out",locVol,color);
	
	inv_param.verbosity=get_quda_verbosity();
	setVerbosityQuda(QUDA_SUMMARIZE,"# QUDA: ",stdout);
	
	/////////////////////////// gauge params ////////////////////////////////
	
	gauge_param=newQudaGaugeParam();
	
	for(int mu=0;mu<NDIM;mu++)
	  gauge_param.X[mu]=locSize[std::array<int,NDIM>{1,2,3,0}[mu]];
	
	gauge_param.anisotropy=1.0;
	
	gauge_param.type=QUDA_WILSON_LINKS;
	gauge_param.gauge_order=QUDA_QDP_GAUGE_ORDER;
	
	gauge_param.t_boundary=QUDA_ANTI_PERIODIC_T;
	
	gauge_param.cpu_prec=QUDA_DOUBLE_PRECISION;
	gauge_param.cuda_prec=QUDA_DOUBLE_PRECISION;
	gauge_param.cuda_prec_sloppy=QUDA_SINGLE_PRECISION;
	gauge_param.cuda_prec_precondition=QUDA_HALF_PRECISION; //check
	gauge_param.cuda_prec_refinement_sloppy=QUDA_SINGLE_PRECISION;
	
	gauge_param.reconstruct=QUDA_RECONSTRUCT_12;
	gauge_param.reconstruct_sloppy=QUDA_RECONSTRUCT_8;
	gauge_param.reconstruct_precondition=QUDA_RECONSTRUCT_8;
	gauge_param.reconstruct_refinement_sloppy=QUDA_RECONSTRUCT_8;
	gauge_param.staggered_phase_type=QUDA_STAGGERED_PHASE_NO;
	gauge_param.staggered_phase_applied=false;//true;
	
	gauge_param.gauge_fix=QUDA_GAUGE_FIXED_NO;
	
	gauge_param.ga_pad=0;
	for(int mu=0;mu<NDIM;mu++)
	  {
	    int surf_size=locVol/locSize[mu]/2;
	    gauge_param.ga_pad=std::max(gauge_param.ga_pad,surf_size);
	  }
	
	inv_param=newQudaInvertParam();
	
	inv_mg_param=newQudaInvertParam();
	quda_mg_param=newQudaMultigridParam();
	quda_mg_param.invert_param=&inv_mg_param;
	
	for(int level=0;level<QUDA_MAX_MG_LEVEL;level++)
	  mg_eig_param[level]=newQudaEigParam();
	
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
	if(quda_mg_preconditioner!=nullptr)
	  {
	    destroyMultigridQuda(quda_mg_preconditioner);
	    quda_mg_preconditioner=nullptr;
	  }
	
	for(int mu=0;mu<NDIM;mu++)
	  nissa_free(quda_conf[mu]);
	
	nissa_free(spincolor_in);
	nissa_free(spincolor_out);
	
	nissa_free(color_in);
	nissa_free(color_out);
	
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
    const double load_clover_time=take_time();
    freeCloverQuda();
    loadCloverQuda(nullptr,nullptr,inv_param);
    master_printf("Time for loadCloverQuda: %lg\n",take_time()-load_clover_time);
  }
  
  /// Reorder conf into QUDA format
  void remap_nissa_to_quda(quda_conf_t out,quad_su3 *in)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	const int iquda=quda_of_loclx[ivol];
	
	for(int nu=0;nu<NDIM;nu++)
	  {
	    const int out_dir=(nu+NDIM-1)%NDIM;
	    su3_copy(out[out_dir][iquda],in[ivol][nu]);
	  }
      }
    NISSA_PARALLEL_LOOP_END;
    
    THREAD_BARRIER();
  }
  
  /// Reorder conf into QUDA format
  void remap_nissa_to_quda(quda_conf_t out,eo_ptr<quad_su3> in)
  {
    for(int par=0;par<2;par++)
      NISSA_PARALLEL_LOOP(ivolh,0,locVolh)
	{
	  const int ivol=loclx_of_loceo[par][ivolh];
	  const int iquda=quda_of_loclx[ivol];
	
	for(int nu=0;nu<NDIM;nu++)
	  {
	    const int out_dir=(nu+NDIM-1)%NDIM;
	    su3_copy(out[out_dir][iquda],in[par][ivolh][nu]);
	  }
      }
    NISSA_PARALLEL_LOOP_END;
    
    THREAD_BARRIER();
  }
  
  /// Reorder spincolor to QUDA format
  void remap_nissa_to_quda(spincolor *out,spincolor *in)
  {
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	const int iquda=quda_of_loclx[ivol];
	spincolor_copy(out[iquda],in[ivol]);
      }
    NISSA_PARALLEL_LOOP_END;
    
    THREAD_BARRIER();
  }
  
  /// Reorder color to QUDA format
  void remap_nissa_to_quda(color *out,eo_ptr<color> in)
  {
    for(int par=0;par<2;par++)
      NISSA_PARALLEL_LOOP(ivolh,0,locVolh)
	{
	  const int ivol=loclx_of_loceo[par][ivolh];
	  const int iquda=quda_of_loclx[ivol];
	  color_copy(out[iquda],in[par][ivolh]);
      }
    NISSA_PARALLEL_LOOP_END;
    
    THREAD_BARRIER();
  }
  
  /// Reorder spincolor from QUDA format
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
  
  /// Reorder color from QUDA format
  void remap_quda_to_nissa(eo_ptr<color> out,color *in)
  {
    NISSA_PARALLEL_LOOP(iquda,0,locVol)
      {
	const int ivol=loclx_of_quda[iquda];
	const int par=loclx_parity[ivol];
	const int ivolh=loceo_of_loclx[ivol];
	color_copy(out[par][ivolh],in[iquda]);
	
#ifdef DEBUG_QUDA
	for(int ic=0;ic<NCOL;ic++)
	  for(int ri=0;ri<2;ri++)
	    {
	      const double& f=in[iquda][ic][ri];
	      if(fabs(f))
		printf("REMA %d %d %d %d %d %lg\n",iquda,par,ivol,ic,ri,f);
	    }
#endif
      }
    NISSA_PARALLEL_LOOP_END;
    
    for(int par=0;par<2;par++)
      set_borders_invalid(out[par]);
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
  
  void set_base_inverter_pars()
  {
    inv_param.dagger=QUDA_DAG_NO;
    inv_param.mass_normalization=QUDA_MASS_NORMALIZATION;
    inv_param.solver_normalization=QUDA_DEFAULT_NORMALIZATION;
    
    inv_param.cpu_prec=QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec=QUDA_DOUBLE_PRECISION;
    inv_param.cuda_prec_sloppy=QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_refinement_sloppy=QUDA_SINGLE_PRECISION;
    inv_param.cuda_prec_precondition=QUDA_HALF_PRECISION;
    inv_param.cuda_prec_eigensolver=QUDA_SINGLE_PRECISION;
    
    inv_param.clover_cpu_prec=QUDA_DOUBLE_PRECISION;
    inv_param.clover_cuda_prec=QUDA_DOUBLE_PRECISION;
    inv_param.clover_cuda_prec_sloppy=QUDA_SINGLE_PRECISION;
    inv_param.clover_cuda_prec_refinement_sloppy=QUDA_SINGLE_PRECISION;
    inv_param.clover_cuda_prec_precondition=QUDA_HALF_PRECISION;
    inv_param.clover_cuda_prec_eigensolver=QUDA_SINGLE_PRECISION;
    
    inv_param.preserve_source=QUDA_PRESERVE_SOURCE_NO;
    inv_param.dirac_order=QUDA_DIRAC_ORDER;
    
    inv_param.input_location=QUDA_CPU_FIELD_LOCATION;
    inv_param.output_location=QUDA_CPU_FIELD_LOCATION;
    
    //inv_param.tune=QUDA_TUNE_YES;
    
    inv_param.sp_pad=0;
    inv_param.cl_pad=0;
    
    inv_param.Ls=1;
    
    inv_param.verbosity=get_quda_verbosity();
    
    inv_param.residual_type=QUDA_L2_RELATIVE_RESIDUAL;
    inv_param.tol_hq=0.1;
    inv_param.reliable_delta=1e-4;
    inv_param.use_sloppy_partial_accumulator=0;
  }
  
  /// Apply the dirac operator
  void apply_tmD(spincolor *out,quad_su3 *conf,double kappa,double csw,double mu,spincolor *in)
  {
    export_gauge_conf_to_external_lib(conf);
    
#ifdef DEBUG_QUDA
    master_printf("setting pars\n");
#endif
    
    inv_param.kappa=kappa;
    set_base_inverter_pars();
    inv_param.gamma_basis=QUDA_CHIRAL_GAMMA_BASIS;
    if(csw)
      {
	inv_param.dslash_type=QUDA_TWISTED_CLOVER_DSLASH;
	inv_param.clover_order=QUDA_PACKED_CLOVER_ORDER;
	inv_param.clover_coeff=csw*kappa;
	// inv_param.clover_cpu_prec=QUDA_DOUBLE_PRECISION;
	// inv_param.clover_cuda_prec=QUDA_DOUBLE_PRECISION;
	// inv_param.cl_pad=0;
	// inv_param.dirac_order=QUDA_DIRAC_ORDER;
	
	loadCloverQuda(NULL,NULL,&inv_param);
      }
    else
      {
	inv_param.dslash_type=QUDA_TWISTED_MASS_DSLASH;
      }
    
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
  
  void set_inverter_pars(const double& kappa,const double& csw,const double& mu,const int& niter,const double& residue,const bool& exported)
  {
    inv_param.kappa=kappa;
    
    if(csw>0.0)
      {
	inv_param.dslash_type=QUDA_TWISTED_CLOVER_DSLASH;
	inv_param.matpc_type=QUDA_MATPC_EVEN_EVEN;
	inv_param.clover_order=QUDA_PACKED_CLOVER_ORDER;
	inv_param.clover_coeff=csw*kappa;
	inv_param.compute_clover=0;
	inv_param.compute_clover_inverse=0;
      }
    else
      {
	inv_param.dslash_type=QUDA_TWISTED_MASS_DSLASH;
	inv_param.matpc_type=QUDA_MATPC_EVEN_EVEN;//QUDA_MATPC_EVEN_EVEN_ASYMMETRIC; //maybe even even
	inv_param.clover_coeff=0;
	inv_param.compute_clover=0;
	inv_param.compute_clover_inverse=0;
      }
    
    inv_param.gamma_basis=QUDA_CHIRAL_GAMMA_BASIS;
    inv_param.solution_type=QUDA_MAT_SOLUTION;
    
    inv_param.inv_type=QUDA_CG_INVERTER;
    inv_param.solve_type=QUDA_NORMERR_PC_SOLVE;
    
    //minus due to different gamma5 definition
    inv_param.mu=-mu;
    inv_param.epsilon=0.0;
    
    inv_param.twist_flavor=QUDA_TWIST_SINGLET;
    inv_param.tol=sqrt(residue);
    inv_param.maxiter=niter;
    inv_param.pipeline=0;
    inv_param.gcrNkrylov=24;
    
    // domain decomposition preconditioner parameters
    inv_param.inv_type_precondition=QUDA_CG_INVERTER;
    //inv_param.schwarz_type=QUDA_ADDITIVE_SCHWARZ;
    inv_param.precondition_cycle=1;
    inv_param.tol_precondition=0.1;
    inv_param.maxiter_precondition=10;
    inv_param.verbosity_precondition=get_quda_verbosity();
    
    inv_param.omega=1.0;
    
    if(exported and
       csw)
      load_clover_term(&inv_param);
    
    if(multiGrid::use_multiGrid)
      {
	inv_param.verbosity=QUDA_SUMMARIZE;//VERBOSE;
	
	// coarsening does not support QUDA_MATPC_EVEN_EVEN_ASYMMETRIC
	if(inv_param.matpc_type==QUDA_MATPC_EVEN_EVEN_ASYMMETRIC)
	  inv_param.matpc_type=QUDA_MATPC_EVEN_EVEN;
	
	inv_param.inv_type=QUDA_GCR_INVERTER;
	inv_param.gcrNkrylov=10; //from default in read_input.l
	inv_param.inv_type_precondition=QUDA_MG_INVERTER;
	inv_param.schwarz_type=QUDA_ADDITIVE_SCHWARZ;
	inv_param.reliable_delta=1e-10;
	inv_param.reliable_delta_refinement=1e-10;
	inv_param.precondition_cycle=1;
	inv_param.tol_precondition=1e-1;
	inv_param.maxiter_precondition=1;
	inv_param.gamma_basis=QUDA_CHIRAL_GAMMA_BASIS;
	inv_param.solve_type=QUDA_DIRECT_SOLVE;
	
	inv_param.omega=1.0;
	
	inv_mg_param=inv_param;
	inv_mg_param.preconditioner=nullptr;
	inv_mg_param.inv_type=QUDA_GCR_INVERTER;
	inv_mg_param.inv_type_precondition=QUDA_INVALID_INVERTER;
	inv_mg_param.maxiter=1000;
	inv_mg_param.solve_type=QUDA_DIRECT_SOLVE;
	inv_mg_param.verbosity=QUDA_SUMMARIZE;//VERBOSE;
	inv_mg_param.residual_type=QUDA_L2_RELATIVE_RESIDUAL;
	inv_mg_param.preserve_source=QUDA_PRESERVE_SOURCE_NO;
	inv_mg_param.gamma_basis=QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
	inv_mg_param.dirac_order=QUDA_DIRAC_ORDER;
	inv_mg_param.input_location=QUDA_CPU_FIELD_LOCATION;
	inv_mg_param.output_location=QUDA_CPU_FIELD_LOCATION;
	inv_mg_param.solution_type=QUDA_MAT_SOLUTION;
	inv_mg_param.dagger=QUDA_DAG_NO;
	
	quda_mg_param.setup_type=QUDA_NULL_VECTOR_SETUP;
	quda_mg_param.pre_orthonormalize=QUDA_BOOLEAN_NO;
	quda_mg_param.post_orthonormalize=QUDA_BOOLEAN_YES;
	
	const int& nlevels=multiGrid::nlevels;
	quda_mg_param.n_level=nlevels;
	
	for(int level=0;level<nlevels;level++)
	  {
	    // set file i/o parameters
	    strcpy(quda_mg_param.vec_infile[level],"");
	    strcpy(quda_mg_param.vec_outfile[level],"");
	    
	    for(int mu=0;mu<NDIM;mu++)
	      {
		const int ext_block_size_mu=multiGrid::block_size[level][mu];
		int& geo_block_size_mu=quda_mg_param.geo_block_size[level][mu];
		
		const int nu=
		  (mu+1)%NDIM;
		
		int extent=glbSize[nu];
		
		//determine how many lattice sites remain at the current level
		for(int k=level;k>0;k--)
		  extent=extent/quda_mg_param.geo_block_size[k-1][mu];
		
		//for the coarsest level, the block size is always set to 1
		if(level==nlevels-1)
		  quda_mg_param.geo_block_size[level][mu]=1;
		else
		  {
		    //the block size for this level and dimension has been set non-zero in the input file
		    //we respect this no matter what
		    if(ext_block_size_mu!=0)
		      quda_mg_param.geo_block_size[level][mu]=ext_block_size_mu;
		    else
		      {
			// on all levels, we try to use a block size of 4^4 and compute the
			// number of fine or aggregate lattice sites on a given level,
			// resulting in block sizes of:
			// - 4 if the extent is larger or equal to 16 and
			// - 2 otherwise
			// When an extent is divisible by three, smaller or equal to 24 and when we're
			// not on the finest grid [and the user has explicitly enabled support
			// for these block lengths  (and therefore also adjusted QUDA to instantiate them)],
			// we use a block length of 3.
			// If aggregation using an even number of lattice points (if size 3 is disabled)
			// is not possible or if the extent is 1 or some divisible only by some prime number
			// other than 3 or 2, we use a block size of 1
			int even_block_size=4;
			if(extent<16)
			  even_block_size=2;
			
			// special treatment of size 24 lattice extents on the fine grid
			if(extent<=24 and extent%3==0)// and quda_input.mg_enable_size_three_blocks)
			  geo_block_size_mu=3;
			else
			  if(extent%even_block_size==0)
			    geo_block_size_mu=even_block_size;
			  else
			    geo_block_size_mu=1;
		      }
		  }
		
		verbosity_lv1_master_printf("# QUDA: MG level %d, extent of (xyzt) dim %d: %d\n",level,mu,extent);
		verbosity_lv1_master_printf("# QUDA: MG aggregation size set to: %d\n",geo_block_size_mu);
		
		//check that all lattice extents are even after blocking on all levels
		if(level < nlevels-1 and (extent/geo_block_size_mu%2))
		  crash("MG level %d, dim %d (xyzt) has extent %d. Block size of %d would result "
			"in odd extent on level %d, aborting!\n"
			"Adjust your block sizes or parallelisation, all local lattice extents on all levels must be even!\n",
			level,mu,extent, geo_block_size_mu,level+1);
		
	      }
	    
	    /// Default value from: https://github.com/lattice/quda/wiki/Multigrid-Solver
	    
	    quda_mg_param.verbosity[level]=get_quda_verbosity();
	    quda_mg_param.precision_null[level]=QUDA_HALF_PRECISION;
	    quda_mg_param.setup_inv_type[level]=QUDA_CG_INVERTER;//QUDA_BICGSTAB_INVERTER or QUDA_CG_INVERTER generally preferred
	    
	    quda_mg_param.num_setup_iter[level]=1;  //Experimental - keep this set to 1
	    quda_mg_param.setup_tol[level]=5e-7;    //Usually around 5e-6 is good.
	    quda_mg_param.setup_maxiter[level]=1000; //500-1000 should work for most systems
	    // If doing twisted mass, we can scale the twisted mass on the coarser grids
	    // which significantly increases speed of convergence as a result of making
	    // the coarsest grid solve a lot better conditioned.
	    // Dean Howarth has some RG arguments on why the coarse mass parameter should be
	    // rescaled for the coarse operator to be optimal.
	    if(fabs(inv_param.mu)>0)
	      {
		quda_mg_param.mu_factor[level]=multiGrid::mu_factor[level];
		master_printf("# QUDA: MG setting coarse mu scaling factor on level %d to %lf\n", level, quda_mg_param.mu_factor[level]);
	      }
	    
	    //Set for all levels except 0. Suggest using QUDA_GCR_INVERTER on all intermediate grids and QUDA_CA_GCR_INVERTER on the bottom.
	    quda_mg_param.coarse_solver[level]=(level+1==nlevels)?QUDA_CA_GCR_INVERTER:QUDA_GCR_INVERTER;
	    constexpr double t[3]={0.25,0.22,0.46};
	    quda_mg_param.coarse_solver_tol[level]=t[level];          //Suggest setting each level to 0.25
	    quda_mg_param.coarse_solver_maxiter[level]=100;//(level+1==nlevels)?50:100;        //Suggest setting in the range 8-100
	    quda_mg_param.spin_block_size[level]=(level==0)?2:1;  //2 for level 0, and 1 thereafter
	    quda_mg_param.n_vec[level]=(level==1)?32:24;          //24 or 32 is supported presently
	    quda_mg_param.nu_pre[level]=nissa::multiGrid::nu_pre[level];            //Suggest setting to 0
	    quda_mg_param.nu_post[level]=nissa::multiGrid::nu_post[level];          //Suggest setting to 8
	    
	    //Always set to QUDA_MG_CYCLE_RECURSIVE (this sets the MG cycles to be a K-cycle which is generally superior to a V-cycle for non-Hermitian systems)
	    quda_mg_param.cycle_type[level]=QUDA_MG_CYCLE_RECURSIVE;
	    //Set to QUDA_CUDA_FIELD_LOCATION for all levels
	    quda_mg_param.location[level]=QUDA_CUDA_FIELD_LOCATION;
	    quda_mg_param.setup_location[level]=QUDA_CUDA_FIELD_LOCATION;
	    
	    quda_mg_param.preserve_deflation=QUDA_BOOLEAN_FALSE;
	    quda_mg_param.smoother[level]=(level+1==nlevels)?QUDA_MR_INVERTER:QUDA_CA_GCR_INVERTER;
	    quda_mg_param.smoother_tol[level]=(level+1==nlevels)?0.25:0.22;                 //Suggest setting each level to 0.25
	    quda_mg_param.smoother_schwarz_cycle[level]=1;          //Experimental, set to 1 for each level
	    //Suggest setting to QUDA_DIRECT_PC_SOLVE for all levels
	    quda_mg_param.smoother_solve_type[level]=QUDA_DIRECT_PC_SOLVE;
	    //Experimental, set to QUDA_INVALID_SCHWARZ for each level unless you know what you're doing
	    quda_mg_param.smoother_schwarz_type[level]=QUDA_INVALID_SCHWARZ;
	    //quda_mg_param.smoother_halo_precision[level]=QUDA_HALF_PRECISION;
	    
	    // when the Schwarz-alternating smoother is used, this can be set to NO, otherwise it must be YES
	    quda_mg_param.global_reduction[level]=QUDA_BOOLEAN_YES;
	    
	    // set to QUDA_MAT_SOLUTION to inject a full field into coarse grid
	    // set to QUDA_MATPC_SOLUTION to inject single parity field into coarse grid
	    // if we are using an outer even-odd preconditioned solve, then we
	    // use single parity injection into the coarse grid
	    quda_mg_param.coarse_grid_solution_type[level]=
	      QUDA_MATPC_SOLUTION;
	    //(inv_param.solve_type==QUDA_DIRECT_PC_SOLVE?QUDA_MATPC_SOLUTION:QUDA_MAT_SOLUTION);
	    quda_mg_param.omega[level]=0.85;  //Set to 0.8-1.0
	    
	    quda_mg_param.location[level]=QUDA_CUDA_FIELD_LOCATION;
	    
	    quda_mg_param.setup_ca_basis[level]     =QUDA_POWER_BASIS;
	    quda_mg_param.setup_ca_basis_size[level]=4;//(level+1==nlevels)?10:4;
	    quda_mg_param.setup_ca_lambda_min[level]=0.0;
	    quda_mg_param.setup_ca_lambda_max[level]=-1.0;
	    
	    quda_mg_param.coarse_solver_ca_basis[level]=QUDA_POWER_BASIS;
	    quda_mg_param.coarse_solver_ca_basis_size[level]=(level+1==nlevels)?10:4;
	    quda_mg_param.coarse_solver_ca_lambda_min[level]=0.0;
	    quda_mg_param.coarse_solver_ca_lambda_max[level]=-1.0;
	    
	    // set the MG EigSolver parameters, almost equivalent to
	    // setEigParam from QUDA's multigrid_invert_test, except
	    // for cuda_prec_ritz (on 20190822)
	    if(level+1==nlevels and multiGrid::use_deflated_solver)
	      {
		quda_mg_param.use_eig_solver[level]=QUDA_BOOLEAN_YES;
		mg_eig_param[level].eig_type=QUDA_EIG_TR_LANCZOS;
		mg_eig_param[level].spectrum=QUDA_SPECTRUM_SR_EIG;
		
		if((mg_eig_param[level].eig_type==QUDA_EIG_TR_LANCZOS or
		    mg_eig_param[level].eig_type==QUDA_EIG_IR_ARNOLDI)
		   and not(mg_eig_param[level].spectrum==QUDA_SPECTRUM_LR_EIG or
			   mg_eig_param[level].spectrum==QUDA_SPECTRUM_SR_EIG))
		  crash("ERROR: MG level %d: Only real spectrum type (LR or SR)"
			"can be passed to the a Lanczos type solver!\n",
			level);
		
		using nissa::multiGrid::nEigenvectors;
		
		mg_eig_param[level].n_ev=nEigenvectors;
		mg_eig_param[level].n_kr=nEigenvectors*1.5;
		mg_eig_param[level].n_conv=nEigenvectors;
		mg_eig_param[level].require_convergence=QUDA_BOOLEAN_TRUE;
		
		mg_eig_param[level].tol=1e-4;
		mg_eig_param[level].check_interval=5;
		mg_eig_param[level].max_restarts=10;
		mg_eig_param[level].cuda_prec_ritz=QUDA_DOUBLE_PRECISION;
		
		mg_eig_param[level].compute_svd=QUDA_BOOLEAN_FALSE;
		mg_eig_param[level].use_norm_op=QUDA_BOOLEAN_TRUE;
		mg_eig_param[level].use_dagger=QUDA_BOOLEAN_FALSE;
		mg_eig_param[level].use_poly_acc=QUDA_BOOLEAN_TRUE;
		mg_eig_param[level].poly_deg=100;
		mg_eig_param[level].a_min=6e-2;
		mg_eig_param[level].a_max=8.0;
		
		// set file i/o parameters
		// Give empty strings, Multigrid will handle IO.
		strcpy(mg_eig_param[level].vec_infile, "");
		strcpy(mg_eig_param[level].vec_outfile, "");
		strncpy(mg_eig_param[level].QUDA_logfile, "quda_eig.log", 512);
		
		quda_mg_param.eig_param[level]=&(mg_eig_param[level]);
	      }
	    else
	      {
		quda_mg_param.eig_param[level]=nullptr;
		quda_mg_param.use_eig_solver[level]=QUDA_BOOLEAN_NO;
	      }
	  }
	
	quda_mg_param.compute_null_vector=QUDA_COMPUTE_NULL_VECTOR_YES;
	quda_mg_param.generate_all_levels=QUDA_BOOLEAN_YES;
	
	// quda_mg_param.run_low_mode_check=QUDA_BOOLEAN_TRUE;//quda_input.mg_run_low_mode_check;
	// quda_mg_param.run_oblique_proj_check=QUDA_BOOLEAN_TRUE;
	// quda_mg_param.run_verify=QUDA_BOOLEAN_TRUE;
	quda_mg_param.run_low_mode_check=QUDA_BOOLEAN_FALSE;//quda_input.mg_run_low_mode_check;
	quda_mg_param.run_oblique_proj_check=QUDA_BOOLEAN_FALSE;
	quda_mg_param.run_verify=QUDA_BOOLEAN_FALSE;
	quda_mg_param.preserve_deflation=QUDA_BOOLEAN_FALSE;
      }
  }
  
  void setup_quda_multigrid()
  {
    bool& setup_valid=multiGrid::setup_valid;
    if(not setup_valid)
      {
	master_printf("QUDA multigrid setup not valid\n");
	
	if(quda_mg_preconditioner!=nullptr)
	  destroyMultigridQuda(quda_mg_preconditioner);
	
	master_printf("mg setup due:\n");
	quda_mg_preconditioner=newMultigridQuda(&quda_mg_param);
	master_printf("mg setup done!\n");
	inv_param.preconditioner=quda_mg_preconditioner;
	
	setup_valid=true;
      }
    // else
    //   {
    // 	master_printf("Updating mg\n");
    // 	updateMultigridQuda(quda_mg_preconditioner,&quda_mg_param);
    //   }
  }
  
  void sanfoPrint(QudaMultigridParam& i)
  {
    int nlev=i.n_level;
    printf("n_level: %d\n",i.n_level);
    for(int ilev=0;ilev<nlev;ilev++)
      for(int idim=0;idim<4;idim++)
	printf("geo_block_size: %d\n",i.geo_block_size[ilev][idim]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("spin_block_size: %d\n",i.spin_block_size[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("n_vec: %d\n",i.n_vec[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("precision_null: %d\n",i.precision_null[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("n_block_ortho: %d\n",i.n_block_ortho[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("verbosity: %d\n",i.verbosity[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("setup_inv_type: %d\n",i.setup_inv_type[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("num_setup_iter: %d\n",i.num_setup_iter[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("setup_tol: %lg\n",i.setup_tol[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("setup_maxiter: %d\n",i.setup_maxiter[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("setup_maxiter_refresh: %d\n",i.setup_maxiter_refresh[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("setup_ca_basis: %d\n",i.setup_ca_basis[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("setup_ca_basis_size: %d\n",i.setup_ca_basis_size[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("setup_ca_lambda_min: %lg\n",i.setup_ca_lambda_min[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("setup_ca_lambda_max: %lg\n",i.setup_ca_lambda_max[ilev]);
    printf("setup_type: %d\n",i.setup_type);
    printf("pre_orthonormalize: %d\n",i.pre_orthonormalize);
    printf("post_orthonormalize: %d\n",i.post_orthonormalize);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("coarse_solver: %d\n",i.coarse_solver[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("coarse_solver_tol: %lg\n",i.coarse_solver_tol[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("coarse_solver_maxiter: %d\n",i.coarse_solver_maxiter[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("coarse_solver_ca_basis: %d\n",i.coarse_solver_ca_basis[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("coarse_solver_ca_basis_size: %d\n",i.coarse_solver_ca_basis_size[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("coarse_solver_ca_lambda_min: %lg\n",i.coarse_solver_ca_lambda_min[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("coarse_solver_ca_lambda_max: %lg\n",i.coarse_solver_ca_lambda_max[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("smoother: %d\n",i.smoother[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("smoother_tol: %lg\n",i.smoother_tol[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("nu_pre: %d\n",i.nu_pre[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("nu_post: %d\n",i.nu_post[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("omega: %lg\n",i.omega[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("smoother_halo_precision: %d\n",i.smoother_halo_precision[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("smoother_schwarz_type: %d\n",i.smoother_schwarz_type[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("smoother_schwarz_cycle: %d\n",i.smoother_schwarz_cycle[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("coarse_grid_solution_type: %d\n",i.coarse_grid_solution_type[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("smoother_solve_type: %d\n",i.smoother_solve_type[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("cycle_type: %d\n",i.cycle_type[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("global_reduction: %d\n",i.global_reduction[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("location: %d\n",i.location[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("setup_location: %d\n",i.setup_location[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("use_eig_solver: %d\n",i.use_eig_solver[ilev]);
    printf("setup_minimize_memory: %d\n",i.setup_minimize_memory);
    printf("compute_null_vector: %d\n",i.compute_null_vector);
    printf("generate_all_levels: %d\n",i.generate_all_levels);
    printf("run_verify: %d\n",i.run_verify);
    printf("run_low_mode_check: %d\n",i.run_low_mode_check);
    printf("run_oblique_proj_check: %d\n",i.run_oblique_proj_check);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("vec_load: %d\n",i.vec_load[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("vec_infile: %s\n",i.vec_infile[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("vec_store: %d\n",i.vec_store[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("vec_outfile: %s\n",i.vec_outfile[ilev]);
    printf("coarse_guess: %d\n",i.coarse_guess);
    printf("preserve_deflation: %d\n",i.preserve_deflation);
    printf("gflops: %lg\n",i.gflops);
    printf("secs: %lg\n",i.secs);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("mu_factor: %lg\n",i.mu_factor[ilev]);
    for(int ilev=0;ilev<nlev;ilev++)
      printf("transfer_type: %d\n",i.transfer_type[ilev]);
    printf("use_mma: %d\n",i.use_mma);
    printf("thin_update_only: %d\n",i.thin_update_only);
  }

  void sanfoPrint(QudaInvertParam& i)
  {
    printf("input_location: %d\n",i.input_location);
    printf("output_location: %d\n",i.output_location);
    printf("dslash_type: %d\n",i.dslash_type);
    printf("inv_type: %d\n",i.inv_type);
    printf("mass: %lg\n",i.mass);
    printf("kappa: %lg\n",i.kappa);
    printf("m5: %lg\n",i.m5);
    printf("Ls: %d\n",i.Ls);
    printf("eofa_shift: %lg\n",i.eofa_shift);
    printf("eofa_pm: %d\n",i.eofa_pm);
    printf("mq1: %lg\n",i.mq1);
    printf("mq2: %lg\n",i.mq2);
    printf("mq3: %lg\n",i.mq3);
    printf("mu: %lg\n",i.mu);
    printf("epsilon: %lg\n",i.epsilon);
    printf("twist_flavor: %d\n",i.twist_flavor);
    printf("laplace3D: %d\n",i.laplace3D);
    printf("tol: %lg\n",i.tol);
    printf("tol_restart: %lg\n",i.tol_restart);
    printf("tol_hq: %lg\n",i.tol_hq);
    printf("compute_true_res: %d\n",i.compute_true_res);
    printf("true_res: %lg\n",i.true_res);
    printf("true_res_hq: %lg\n",i.true_res_hq);
    printf("maxiter: %d\n",i.maxiter);
    printf("reliable_delta: %lg\n",i.reliable_delta);
    printf("reliable_delta_refinement: %lg\n",i.reliable_delta_refinement);
    printf("use_alternative_reliable: %d\n",i.use_alternative_reliable);
    printf("use_sloppy_partial_accumulator: %d\n",i.use_sloppy_partial_accumulator);
    printf("solution_accumulator_pipeline: %d\n",i.solution_accumulator_pipeline);
    printf("max_res_increase: %d\n",i.max_res_increase);
    printf("max_res_increase_total: %d\n",i.max_res_increase_total);
    printf("max_hq_res_increase: %d\n",i.max_hq_res_increase);
    printf("max_hq_res_restart_total: %d\n",i.max_hq_res_restart_total);
    printf("heavy_quark_check: %d\n",i.heavy_quark_check);
    printf("pipeline: %d\n",i.pipeline);
    printf("num_offset: %d\n",i.num_offset);
    printf("num_src: %d\n",i.num_src);
    printf("num_src_per_sub_partition: %d\n",i.num_src_per_sub_partition);
    for(int idim=0;idim<4;idim++)
      printf("split_grid: %d\n",i.split_grid[idim]);
    printf("overlap: %d\n",i.overlap);
    printf("compute_action: %d\n",i.compute_action);
    printf("solution_type: %d\n",i.solution_type);
    printf("solve_type: %d\n",i.solve_type);
    printf("matpc_type: %d\n",i.matpc_type);
    printf("dagger: %d\n",i.dagger);
    printf("mass_normalization: %d\n",i.mass_normalization);
    printf("solver_normalization: %d\n",i.solver_normalization);
    printf("preserve_source: %d\n",i.preserve_source);
    printf("cpu_prec: %d\n",i.cpu_prec);
    printf("cuda_prec: %d\n",i.cuda_prec);
    printf("cuda_prec_sloppy: %d\n",i.cuda_prec_sloppy);
    printf("cuda_prec_refinement_sloppy: %d\n",i.cuda_prec_refinement_sloppy);
    printf("cuda_prec_precondition: %d\n",i.cuda_prec_precondition);
    printf("cuda_prec_eigensolver: %d\n",i.cuda_prec_eigensolver);
    printf("dirac_order: %d\n",i.dirac_order);
    printf("gamma_basis: %d\n",i.gamma_basis);
    printf("clover_location: %d\n",i.clover_location);
    printf("clover_cpu_prec: %d\n",i.clover_cpu_prec);
    printf("clover_cuda_prec: %d\n",i.clover_cuda_prec);
    printf("clover_cuda_prec_sloppy: %d\n",i.clover_cuda_prec_sloppy);
    printf("clover_cuda_prec_refinement_sloppy: %d\n",i.clover_cuda_prec_refinement_sloppy);
    printf("clover_cuda_prec_precondition: %d\n",i.clover_cuda_prec_precondition);
    printf("clover_cuda_prec_eigensolver: %d\n",i.clover_cuda_prec_eigensolver);
    printf("clover_order: %d\n",i.clover_order);
    printf("use_init_guess: %d\n",i.use_init_guess);
    printf("clover_coeff: %lg\n",i.clover_coeff);
    printf("clover_rho: %lg\n",i.clover_rho);
    printf("compute_clover_trlog: %d\n",i.compute_clover_trlog);
    printf("compute_clover: %d\n",i.compute_clover);
    printf("compute_clover_inverse: %d\n",i.compute_clover_inverse);
    printf("return_clover: %d\n",i.return_clover);
    printf("return_clover_inverse: %d\n",i.return_clover_inverse);
    printf("verbosity: %d\n",i.verbosity);
    printf("sp_pad: %d\n",i.sp_pad);
    printf("cl_pad: %d\n",i.cl_pad);
    printf("iter: %d\n",i.iter);
    printf("gflops: %lg\n",i.gflops);
    printf("secs: %lg\n",i.secs);
    printf("tune: %d\n",i.tune);
    printf("Nsteps: %d\n",i.Nsteps);
    printf("gcrNkrylov: %d\n",i.gcrNkrylov);
    printf("inv_type_precondition: %d\n",i.inv_type_precondition);
    printf("deflate: %d\n",i.deflate);
    printf("dslash_type_precondition: %d\n",i.dslash_type_precondition);
    printf("verbosity_precondition: %d\n",i.verbosity_precondition);
    printf("tol_precondition: %lg\n",i.tol_precondition);
    printf("maxiter_precondition: %d\n",i.maxiter_precondition);
    printf("omega: %lg\n",i.omega);
    printf("ca_basis: %d\n",i.ca_basis);
    printf("ca_lambda_min: %lg\n",i.ca_lambda_min);
    printf("ca_lambda_max: %lg\n",i.ca_lambda_max);
    printf("precondition_cycle: %d\n",i.precondition_cycle);
    printf("schwarz_type: %d\n",i.schwarz_type);
    printf("residual_type: %d\n",i.residual_type);
    printf("cuda_prec_ritz: %d\n",i.cuda_prec_ritz);
    printf("n_ev: %d\n",i.n_ev);
    printf("max_search_dim: %d\n",i.max_search_dim);
    printf("rhs_idx: %d\n",i.rhs_idx);
    printf("deflation_grid: %d\n",i.deflation_grid);
    printf("eigenval_tol: %lg\n",i.eigenval_tol);
    printf("eigcg_max_restarts: %d\n",i.eigcg_max_restarts);
    printf("max_restart_num: %d\n",i.max_restart_num);
    printf("inc_tol: %lg\n",i.inc_tol);
    printf("make_resident_solution: %d\n",i.make_resident_solution);
    printf("use_resident_solution: %d\n",i.use_resident_solution);
    printf("chrono_make_resident: %d\n",i.chrono_make_resident);
    printf("chrono_replace_last: %d\n",i.chrono_replace_last);
    printf("chrono_use_resident: %d\n",i.chrono_use_resident);
    printf("chrono_max_dim: %d\n",i.chrono_max_dim);
    printf("chrono_index: %d\n",i.chrono_index);
    printf("chrono_precision: %d\n",i.chrono_precision);
    printf("extlib_type: %d\n",i.extlib_type);
    printf("native_blas_lapack: %d\n",i.native_blas_lapack);
  }
  
  bool solve_tmD(spincolor *sol,quad_su3 *conf,const double& kappa,const double& csw,const double& mu,const int& niter,const double& residue,spincolor *source)
  {
    const double export_time=take_time();
    const bool exported=export_gauge_conf_to_external_lib(conf);
    master_printf("time to export to external library: %lg s\n",take_time()-export_time);
    
    set_base_inverter_pars();
    
    set_inverter_pars(kappa,csw,mu,niter,residue,exported);
    
#ifdef DEBUG_QUDA
    if(is_master_rank())
      {
	master_printf("--- gauge pars: ---\n");
	printQudaGaugeParam(&gauge_param);
	master_printf("--- inv pars: ---\n");
	printQudaInvertParam(&inv_param);
	
	if(multiGrid::use_multiGrid)
	  {
	    master_printf("--- inv_mg pars: %p kappa %lg mu %lg mass %lg---\n",&inv_mg_param,inv_mg_param.kappa,inv_mg_param.mu,inv_mg_param.mass);
	    printQudaInvertParam(&inv_mg_param);
	    master_printf("--- multigrid pars: %p internal %p ---\n",&quda_mg_param,quda_mg_param.invert_param);
	    printQudaMultigridParam(&quda_mg_param);
	    master_printf("AAAAAA\n");
	    if(is_master_rank())
	      sanfoPrint(quda_mg_param);
	    master_printf("AAAAAA\n");
	    if(is_master_rank())
	      sanfoPrint(*quda_mg_param.invert_param);
	    master_printf("AAAAAA\n");
	    
	    master_printf("-- -eig pars: ---\n");
	    printQudaEigParam(mg_eig_param+multiGrid::nlevels-1);
	  }
      }
#endif
    
    if(multiGrid::use_multiGrid)
      setup_quda_multigrid();
    
    const double remap_in_time=take_time();
    remap_nissa_to_quda(spincolor_in,source);
    master_printf("time to remap rhs to quda: %lg s\n",take_time()-remap_in_time);
    
    const double solution_time=take_time();
    invertQuda(spincolor_out,spincolor_in,&inv_param);
    master_printf("Solution time: %lg s\n",take_time()-solution_time);
    
    master_printf("# QUDA solved in: %i iter / %g secs=%g Gflops\n",inv_param.iter,inv_param.secs,inv_param.gflops/inv_param.secs);
    
    const double remap_out_time=take_time();
    remap_quda_to_nissa(sol,spincolor_out);
    master_printf("time to remap solution from quda: %lg s\n",take_time()-remap_out_time);
    
    // Might return actual result of the convergence with some proper error handling?
    return true;
  }
  
  bool solve_stD(eo_ptr<color> sol,eo_ptr<quad_su3> conf,const double& mass,const int& niter,const double& residue,eo_ptr<color> source)
  {
    gauge_param.reconstruct=
      gauge_param.reconstruct_sloppy=
      gauge_param.reconstruct_precondition=
      gauge_param.reconstruct_refinement_sloppy=
      QUDA_RECONSTRUCT_NO;
    const double export_time=take_time();
    
    add_or_rem_stagphases_to_conf(conf);
    const bool exported=export_gauge_conf_to_external_lib(conf);
    add_or_rem_stagphases_to_conf(conf);
    master_printf("time to export (%d) to external library: %lg s\n",exported,take_time()-export_time);
    
    set_base_inverter_pars();
    
    inv_param.dslash_type=QUDA_STAGGERED_DSLASH;
    inv_param.gamma_basis=QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
    inv_param.matpc_type=QUDA_MATPC_EVEN_EVEN;
    inv_param.solution_type=QUDA_MAT_SOLUTION;
    
    inv_param.inv_type=QUDA_CG_INVERTER;
    inv_param.solve_type=QUDA_DIRECT_PC_SOLVE;
    inv_param.dagger=QUDA_DAG_NO;
    
    //minus due to different gamma5 definition
    inv_param.mass=mass;
    
    inv_param.tol=sqrt(residue);
    inv_param.maxiter=niter;
    inv_param.Ls=1;
    
    inv_param.verbosity=get_quda_verbosity();
    
    remap_nissa_to_quda(color_in,source);
    
    //invertQuda(color_out,color_in,&inv_param);
    MatQuda(color_out,color_in,&inv_param);
    
    master_printf("# QUDA solved in: %i iter / %g secs=%g Gflops\n",inv_param.iter,inv_param.secs,inv_param.gflops/inv_param.secs);
    
    remap_quda_to_nissa(sol,color_out);
    
    return true;
  }
}
