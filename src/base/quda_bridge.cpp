#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#define EXTERN_QUDA_BRIDGE
# include "quda_bridge.hpp"

#include "base/cuda.hpp"
#include "io/checksum.hpp"
#include "base/export_conf_to_external_solver.hpp"
#include "base/multiGridParams.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "metaprogramming/hasMember.hpp"
#include "new_types/su3_op.hpp"
#include "new_types/custom_real_numb.hpp"
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

#ifdef READY_TO_DEL
namespace quda
{
  /// Make visible quda headers not exported
  void setDiracParam(DiracParam &diracParam,QudaInvertParam *inv_param,const bool pc);
}
#endif

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

namespace nissa
{
  PROVIDE_HAS_MEMBER(cl_pad);
  PROVIDE_HAS_MEMBER(sp_pad);
  PROVIDE_HAS_MEMBER(gflops);
  PROVIDE_HAS_MEMBER(n_vec_batch);
  PROVIDE_HAS_MEMBER(secs);
  PROVIDE_HAS_MEMBER(true_res);
  PROVIDE_HAS_MEMBER(true_res_hq);
  PROVIDE_HAS_MEMBER(tune);
}

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
  CUDA_MANAGED spincolor *spincolor_in=nullptr;
  CUDA_MANAGED spincolor *spincolor_out=nullptr;
  
  /// Return the rank of the given quda coords
  int get_rank_of_quda_coords(const int *coords,void *fdata)
  {
    /// Coordinates
    int c[NDIM];
    for(int mu=0;mu<NDIM;mu++)
      c[mu]=coords[std::array<int,NDIM>{3,0,1,2}[mu]];
    
    /// Result
    int out;
    MPI_Cart_rank(quda_iface::cart_comm,c,&out);
    
    verbosity_lv3_master_printf("rank %d, {%d %d %d %d}->%d\n",rank,c[0],c[1],c[2],c[3],out);
    
    return out;
  }
  
  const char* getPrecTag(const int precision)
    {
      switch (precision)
	{
	case QUDA_DOUBLE_PRECISION:
	  return "QUDA_DOUBLE_PRECISION";
	case QUDA_SINGLE_PRECISION:
	  return "QUDA_SINGLE_PRECISION";
	case QUDA_HALF_PRECISION:
	  return "QUDA_HALF_PRECISION";
	case QUDA_QUARTER_PRECISION:
	  return "QUDA_QUARTER_PRECISION";
	  break;
	default:
	  return "";
	}
    }
  
  namespace internal
  {
    template <QudaPrecision>
    struct _CustomRealOfQudaPrecision;

#define PROVIDE_CUSTOM_REAL_OF_QUDA_PRECISION(QUDA_ENUM,TYPE)	\
    template <>							\
    struct _CustomRealOfQudaPrecision<QUDA_ENUM>		\
    {								\
      using type=TYPE;						\
    }
    
    PROVIDE_CUSTOM_REAL_OF_QUDA_PRECISION(QUDA_DOUBLE_PRECISION,CustomDouble);
    PROVIDE_CUSTOM_REAL_OF_QUDA_PRECISION(QUDA_SINGLE_PRECISION,CustomFloat);
    PROVIDE_CUSTOM_REAL_OF_QUDA_PRECISION(QUDA_HALF_PRECISION,CustomHalf);
    
#undef PROVIDE_CUSTOM_REAL_OF_QUDA_PRECISION
  }
  
  /// Custom type corresponding to quda precision
  template <QudaPrecision Prec>
  using CustomRealOfQudaPrecision=
    typename internal::_CustomRealOfQudaPrecision<Prec>::type;
  
  // /// Gets the i-th entry of an array v, interpeting entries as Prec type 
  // double getFromCustomPrecArray(const void* v,
  // 				const size_t& i,
  // 				const size_t prec)
  // {
  //   if(prec==8)
  //     return (double)((CustomRealOfQudaPrecision<QUDA_DOUBLE_PRECISION>*)v)[i];
  //   else if(prec==4)
  //     return (double)((CustomRealOfQudaPrecision<QUDA_SINGLE_PRECISION>*)v)[i];
  //   else if(prec==QUDA_HALF_PRECISION)
  //     return __half2float(((__half*)v)[i]);
  //   else
  //     crash("Unknown precision %zu",prec);
    
  //   return 0;
  // }
  
#ifdef QUDASETUPSTORE
  void QudaSetup::restoreOrTakeCopyOfB(const bool takeCopy,
				       std::vector<quda::ColorSpinorField*>& Bdev,
				       const size_t lev)
  {
    using namespace nissa::Robbery;
    
    const size_t nB=Bdev.size();
    const int prec=Bdev[0]->Precision();
    const size_t byteSize=nB?(Bdev[0]->Bytes()):0;
    const size_t ghostSize=nB?Bdev[0]->GhostBytes():0;
    master_printf("B size: %zu bytes, precision %d (%s) for each of the %zu vectors, corresponding to %zu complex, ghost size: %zu\n",byteSize,prec,getPrecTag(prec),nB,byteSize/(2*prec),ghostSize);
    
    if(takeCopy)
      {
	B[lev].resize(nB);
	for(size_t iB=0;iB<nB;iB++)
	  B[lev][iB]=nissa_malloc(("Bi"+std::to_string(iB)).c_str(),byteSize,char);
	allocatedMemory+=nB*byteSize;
      }
    else
      if(const size_t nBHost=B[lev].size();nBHost!=nB)
	crash("B size not matching, this is %zu and device setup is %zu",nBHost,nB);
    
    for(size_t iB=0;iB<nB;iB++)
      {
	restoreOrTakeCopyOfData(B[lev][iB],Bdev[iB]->V(),byteSize,takeCopy);
	
	master_printf("B[%zu] vec of lev %zu %s, first entries: ",iB,lev,takeCopy?"stored":"restored");
	for(int i=0;i<2;i++)
	  master_printf("%lg ",getFromCustomPrecArray((B[lev])[iB],i,prec));
	
	master_printf("\n");
      }
  }
  
  void QudaSetup::restoreOrTakeCopyOfEig(const bool takeCopy,
					 quda::Solver* csv)
  {
    using namespace nissa::Robbery;
    
    quda::Solver* nestedSolver=rob<solver>((quda::PreconditionedSolver*)csv);
    master_printf("nested solver: %p\n",nestedSolver);
    
    auto& eVecsDev=rob<evecs>(nestedSolver);
    const size_t nEig=eVecsDev.size();
    master_printf("nEig: %zu\n",nEig);
    const size_t byteSize=nEig?(eVecsDev[0]->Bytes()):0;
    
    if(takeCopy)
      {
	eVecs.resize(nEig);
	for(size_t iEig=0;iEig<nEig;iEig++)
	  eVecs[iEig]=nissa_malloc(("ei"+std::to_string(iEig)).c_str(),byteSize,char);
	allocatedMemory+=byteSize*nEig;
      }
    else
      if(nEig!=eVecs.size())
	crash("eig size not matching, this is %zu and device setup is %zu",eVecs.size(),nEig);
    
    for(size_t iEig=0;iEig<nEig;iEig++)
      restoreOrTakeCopyOfData(eVecs[iEig],
			      eVecsDev[iEig]->V(),
			      byteSize,
			      takeCopy);
    
    auto& eValsDev=rob<evals>(nestedSolver);
    if(takeCopy)
      eVals=eValsDev;
    else
      eValsDev=eVals;
    
    master_printf("eigenvecs %s, first entries of evals: (%lg,%lg) (%lg,%lg), of eVecs of prec %d: ",takeCopy?"stored":"restored",
		  eVals[0].real(),eVals[0].imag(),
		  eVals[1].real(),eVals[1].imag(),
		  eVecsDev[0]->Precision());
    for(int i=0;i<2;i++)
      master_printf("(%lg,%lg) ",
		    getFromCustomPrecArray(eVecs[i],0+2*i,eVecsDev[0]->Precision()),
		    getFromCustomPrecArray(eVecs[i],1+2*i,eVecsDev[0]->Precision()));
    master_printf("\n");
  }
  
  void QudaSetup::restoreOrTakeCopyOfAllY(const bool takeCopy)
  {
    using namespace nissa::Robbery;
    using namespace quda;
    
    MG* cur=static_cast<multigrid_solver*>(quda_mg_preconditioner)->mg;
    int lev=0;
    
    while(lev<multiGrid::nlevels-1)
      {
	DiracCoarse* dc=static_cast<DiracCoarse*>(rob<diracCoarseSmoother>(cur));
	cudaGaugeField* yd=rob<Y_d>(dc);
	cudaGaugeField* yhat_d=rob<Yhat_d>(dc);
	
	if(takeCopy)
	  {
	    Y[lev]=nissa_malloc("Y",yd->Bytes(),char);
	    Yhat[lev]=nissa_malloc("Yhat",yhat_d->Bytes(),char);
	    allocatedMemory+=yd->Bytes()+yhat_d->Bytes();
	  }
	
	restoreOrTakeCopyOfData(Y[lev],yd->Gauge_p(),yd->Bytes(),takeCopy);
	restoreOrTakeCopyOfData(Yhat[lev],yhat_d->Gauge_p(),yhat_d->Bytes(),takeCopy);
	
	master_printf("%s Y and Yhat, lev %d, size %d\n",takeCopy?"stored":"restored",lev,(int)yd->Bytes());
	for(int i=0;i<10;i++)
	  master_printf("y[%zu]: %.16lg\n",i,getFromCustomPrecArray(Y[lev],i,yd->Precision()));
	cur=rob<coarse>(cur);
	
	lev++;
      }
  }
  
  void QudaSetup::restoreOrTakeCopy(const bool takeCopy)
  {
    using namespace nissa::Robbery;
    using namespace quda;
    
    multigrid_solver* mgs=static_cast<multigrid_solver*>(quda_mg_preconditioner);
    MG* cur=mgs->mg;
    master_printf("/////////////////////////////////////////////////////////////////\n");
    master_printf("/////////////////////////// preverify //////////////////////////////////////\n");
    cur->verify();
    master_printf("/////////////////////////////////////////////////////////////////\n");
    int lev=0;
    
    allocatedMemory=0;
    
    if(takeCopy)
      {
	if(not (B.empty() and Y.empty() and Yhat.empty())) crash("setup already in use!");
	B.resize(multiGrid::nlevels);
	for(auto& p : {&Y,&Yhat})
	  p->resize(multiGrid::nlevels);
      }
    else
      if(B.empty() or Y.empty() or Yhat.empty()) crash("setup not in use!");
    
    restoreOrTakeCopyOfAllY(takeCopy);
    
    master_printf("&mgs->B %p , &mgs->mgParam.B %p\n",&mgs->B,&mgs->mgParam->B);
    restoreOrTakeCopyOfB(takeCopy,mgs->mgParam->B,lev);
    
    lev=1;
    while(lev<multiGrid::nlevels)
      {
	MGParam* mgLevParam=rob<param_coarse>(cur);
	std::vector<ColorSpinorField*>& Bdev=mgLevParam->B;
	restoreOrTakeCopyOfB(takeCopy,Bdev,lev);
	
	//Dirac* dc=rob<diracCoarseSmoother>(cur);
	Solver* csv=rob<coarse_solver>(cur);
	master_printf("csv: %p\n",csv);
	if(csv and quda_mg_param.use_eig_solver[lev]==QUDA_BOOLEAN_YES)
	  {
	    master_printf("Going to to the eig part\n");
	    restoreOrTakeCopyOfEig(takeCopy,csv);
	  }
	
	cur=rob<coarse>(cur);
	
	master_printf("Done with lev %d\n",lev);
	
	lev++;
      }
    
    if(takeCopy)
      master_printf("we have used %zu bytes to store the setup\n",allocatedMemory);
    else
      {
	master_printf("Everything recycled in the mg\n");
	updateMultigridQuda(quda_mg_preconditioner,&quda_mg_param);
	multigrid_solver* mgs=static_cast<multigrid_solver*>(quda_mg_preconditioner);
	MG* cur=mgs->mg;
	master_printf("Let us verify\n");
	cur->verify();
	
	// QudaBoolean& p=quda_mg_param.preserve_deflation; //thin_update_only
	// const QudaBoolean oldP=p;
	// p=QUDA_BOOLEAN_TRUE;
	// updateMultigridQuda(quda_mg_preconditioner,&quda_mg_param);
	// p=oldP;
      }
  }
#endif
  
  /// Set the verbosity
  QudaVerbosity get_verbosity_for_quda()
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
	coords_t periods;
	for(int mu=0;mu<NDIM;mu++)
	  periods[mu]=1;
	MPI_Cart_create(MPI_COMM_WORLD,NDIM,&nrank_dir[0],&periods[0],1,&cart_comm);
	
	master_printf("Initializing QUDA\n");
	
	if(QUDA_VERSION_MAJOR==0 and QUDA_VERSION_MINOR<7)
	  crash("minimum QUDA version required is 0.7.0");
	
	/////////////////////////////////////////////////////////////////
	
	quda_of_loclx=nissa_malloc("quda_of_loclx",locVol,int);
	loclx_of_quda=nissa_malloc("loclx_of_quda",locVol,int);
	
	for(int64_t ivol=0;ivol<locVol;ivol++)
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
	      crash("quda %d remapping to ivol %ld not in range [0,%ld]",quda,ivol,locVol);
	    
	    quda_of_loclx[ivol]=quda;
	    loclx_of_quda[quda]=ivol;
	  }
	
	master_printf("allocating quda_conf\n");
	for(int mu=0;mu<NDIM;mu++)
	  quda_conf[mu]=nissa_malloc("quda_conf",locVol,su3);
	
	spincolor_in=nissa_malloc("spincolor_in",locVol,spincolor);
	spincolor_out=nissa_malloc("spincolor_out",locVol,spincolor);
	
	color_in=nissa_malloc("color_in",locVol,color);
	color_out=nissa_malloc("color_out",locVol,color);
	
	setVerbosityQuda(QUDA_VERBOSE,"# QUDA: ",stdout);
	
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
	
	gauge_param.reconstruct=QUDA_RECONSTRUCT_NO;//QUDA_RECONSTRUCT_12;
	gauge_param.reconstruct_sloppy=QUDA_RECONSTRUCT_NO;
	gauge_param.reconstruct_precondition=QUDA_RECONSTRUCT_NO;
	gauge_param.reconstruct_refinement_sloppy=QUDA_RECONSTRUCT_NO;
	gauge_param.staggered_phase_type=QUDA_STAGGERED_PHASE_NO;
	gauge_param.staggered_phase_applied=false;//true;
	
	gauge_param.gauge_fix=QUDA_GAUGE_FIXED_NO;
	
	gauge_param.ga_pad=0;
	for(int mu=0;mu<NDIM;mu++)
	  {
	    int surf_size=locVol/locSize[mu]/2;
	    gauge_param.ga_pad=std::max(gauge_param.ga_pad,surf_size);
	  }
	
	int grid[NDIM];
	for(int mu=0;mu<NDIM;mu++)
	  grid[mu]=nrank_dir[std::array<int,NDIM>{1,2,3,0}[mu]];
	
	initCommsGridQuda(NDIM,grid,get_rank_of_quda_coords,NULL);
	
	// initQuda(iCudaDevice);
	initQuda(-1);
	
	setVerbosityQuda(QUDA_VERBOSE,"# QUDA: ",stdout);
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
	
#ifdef QUDASETUPSTORE
	quda_iface::qudaSetups.clear();
#endif
	
	//free the clover and gauge conf
	freeGaugeQuda();
#ifndef DYNAMIC_CLOVER
	freeCloverQuda();
#endif
	endQuda();
	
	inited=false;
      }
  }
  
  /// Loads the clover term
  void load_clover_term(QudaInvertParam* inv_param)
  {
    const double load_clover_time=
      take_time();
    freeCloverQuda();
    loadCloverQuda(nullptr,nullptr,inv_param);
    master_printf("Time for loadCloverQuda: %lg\n",take_time()-load_clover_time);
  }
  
  /// Reorder conf into QUDA format
  void remap_nissa_to_quda(quda_conf_t& out,
			   const LxField<quad_su3>& in)
  {
    PAR(0,
	locVol,
	CAPTURE(TO_READ(in),
		out),
	ivol,
	{
	  const int iquda=quda_of_loclx[ivol];
	  
#pragma unroll 4
	  for(int nu=0;nu<NDIM;nu++)
	    {
	      const int out_dir=(nu+NDIM-1)%NDIM;
	      su3_copy(out[out_dir][iquda],in[ivol][nu]);
	    }
	});
  }
  
  /// Reorder conf into QUDA format
  void remap_nissa_to_quda(quda_conf_t& out,
			   const EoField<quad_su3>& in)
  {
    for(int par=0;par<2;par++)
      PAR(0,
	  locVolh,
	  CAPTURE(TO_READ(in),
		  out,
		  par),
	  ivolh,
	  {
	    const int ivol=loclx_of_loceo[par][ivolh];
	    const int iquda=quda_of_loclx[ivol];
	    
#pragma unroll 4
	    for(int nu=0;nu<NDIM;nu++)
	      {
		const int out_dir=(nu+NDIM-1)%NDIM;
		su3_copy(out[out_dir][iquda],in[par][ivolh][nu]);
	      }
	  });
  }
  
  /// Reorder spincolor to QUDA format
  void remap_nissa_to_quda(spincolor* out,
			   const LxField<spincolor>& in)
  {
    PAR(0,
	locVol,
	CAPTURE(TO_READ(in),
		out),
	ivol,
	{
	  const int iquda=quda_of_loclx[ivol];
	  spincolor_copy(out[iquda],in[ivol]);
	});
  }
  
  /// Reorder color to QUDA format
  void remap_nissa_to_quda(color *out,
			   const EoField<color>& in)
  {
    for(int par=0;par<2;par++)
      PAR(0,
	  locVolh,
	  CAPTURE(TO_READ(in),
		  out,
		  par),
	  ivolh,
	  {
	    const int ivol=loclx_of_loceo[par][ivolh];
	    const int iquda=quda_of_loclx[ivol];
	    color_copy(out[iquda],in[par][ivolh]);
	  });
  }
  
  /// Reorder spincolor from QUDA format
  void remap_quda_to_nissa(LxField<spincolor>& out,
			   const spincolor *in)
  {
    PAR(0,
	locVol,
	CAPTURE(TO_WRITE(out),
		in),
	iquda,
	{
	  const int ivol=loclx_of_quda[iquda];
	  spincolor_copy(out[ivol],in[iquda]);
	});
  }
  
  /// Reorder color from QUDA format
  void remap_quda_to_nissa(EoField<color>& out,
			   color* in)
  {
    PAR(0,
	locVol,
	CAPTURE(in,
		TO_WRITE(out)),
	iquda,
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
	});
  }
  
  // /// Sets the sloppy precision
  // void set_sloppy_prec(const QudaPrecision sloppy_prec)
  // {
  //   master_printf("Quda sloppy precision: ");
  
  //   switch(sloppy_prec)
  //     {
  //     case QUDA_DOUBLE_PRECISION:
  // 	master_printf("double");
  // 	break;
  //     case QUDA_SINGLE_PRECISION:
  // 	master_printf("single");
  // 	break;
  //     case QUDA_HALF_PRECISION:
  // 	master_printf("half");
  // 	break;
  //     case QUDA_QUARTER_PRECISION:
  // 	master_printf("quarter");
  // 	break;
  //     case QUDA_INVALID_PRECISION:
  // 	crash("invalid precision");
  //     };
  //   master_printf("\n");
  
  // gauge_param.cuda_prec_sloppy=
  //   gauge_param.cuda_prec_refinement_sloppy=
  //   inv_param.cuda_prec_sloppy=
  //   inv_param.clover_cuda_prec_sloppy=
  //   inv_param.clover_cuda_prec_refinement_sloppy=
  //   sloppy_prec;
  // }

#define PROVIDE_MAYBE_SET(X)					\
  template <typename T,						\
	    typename V>						\
  void maybe_set_ ## X(T& i,					\
		       V&& v)					\
  {								\
    if constexpr(hasMember_ ## X<T>)				\
      i.X=v;							\
  }
  
  PROVIDE_MAYBE_SET(cl_pad);
  PROVIDE_MAYBE_SET(sp_pad);
  
#undef PROVIDE_MAYBE_SET
  
  void set_base_inverter_pars()
  {
    inv_param=newQudaInvertParam();
    
    inv_param.verbosity=get_verbosity_for_quda();
    
    inv_mg_param=newQudaInvertParam();
    quda_mg_param=newQudaMultigridParam();
    quda_mg_param.invert_param=&inv_mg_param;
    
    for(int level=0;level<QUDA_MAX_MG_LEVEL;level++)
      mg_eig_param[level]=newQudaEigParam();
    
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
    
    inv_param.input_location=QUDA_CUDA_FIELD_LOCATION;
    inv_param.output_location=QUDA_CUDA_FIELD_LOCATION;
    
    //inv_param.tune=QUDA_TUNE_YES;
    
    maybe_set_sp_pad(inv_param,0);
    maybe_set_cl_pad(inv_param,0);
    
    inv_param.Ls=1;
    
    inv_param.verbosity=get_verbosity_for_quda();
    
    inv_param.residual_type=QUDA_L2_RELATIVE_RESIDUAL;
    inv_param.tol_hq=0.1;
    inv_param.reliable_delta=nissa::multiGrid::reliable_delta;
    inv_param.use_sloppy_partial_accumulator=0;
    inv_param.chrono_precision=QUDA_SINGLE_PRECISION;
  }
  
  /// Apply the dirac operator
  void apply_tmD(spincolor *out,quad_su3 *conf,double kappa,double csw,double mu,spincolor *in)
  {
    crash("reimplement");
    //     export_gauge_conf_to_external_solver(conf);
    
// #ifdef DEBUG_QUDA
//     master_printf("setting pars\n");
// #endif
    
    // set_base_inverter_pars();
    // inv_param.kappa=kappa;
    // inv_param.gamma_basis=QUDA_CHIRAL_GAMMA_BASIS;
    // if(csw)
    //   {
    // 	inv_param.dslash_type=QUDA_TWISTED_CLOVER_DSLASH;
    // 	inv_param.clover_order=QUDA_PACKED_CLOVER_ORDER;
    // 	inv_param.clover_coeff=csw*kappa;
    // 	// inv_param.clover_cpu_prec=QUDA_DOUBLE_PRECISION;
    // 	// inv_param.clover_cuda_prec=QUDA_DOUBLE_PRECISION;
	
    // 	if constexpr(hasMember_cl_pad<decltype(inv_param)>)
    // 	  inv_param.cl_pad=0;
	
// 	loadCloverQuda(NULL,NULL,&inv_param);
//       }
//     else
//       {
// 	inv_param.dslash_type=QUDA_TWISTED_MASS_DSLASH;
//       }
    
//     inv_param.solution_type=QUDA_MAT_SOLUTION;
    
//     //minus due to different gamma5 definition
//     inv_param.mu=-mu;
//     inv_param.epsilon=0.0;
    
//     inv_param.twist_flavor=QUDA_TWIST_SINGLET;
//     inv_param.Ls=1;
    
//     remap_nissa_to_quda(spincolor_in,in);
//     MatQuda(spincolor_out,spincolor_in,&inv_param);
//     remap_quda_to_nissa(out,spincolor_out);
  }
  
  template <typename T>
  void set_n_vec_batch(T& quda_mg_param)
  {
    if constexpr(hasMember_n_vec_batch<T>)
      for(int level=0;level<multiGrid::nlevels;level++)
	quda_mg_param.n_vec_batch[level]=1;
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
    inv_param.gcrNkrylov=nissa::multiGrid::gcrNkrylov;
    
    // domain decomposition preconditioner parameters
    inv_param.inv_type_precondition=QUDA_CG_INVERTER;
    //inv_param.schwarz_type=QUDA_ADDITIVE_SCHWARZ;
    inv_param.precondition_cycle=1;
    inv_param.tol_precondition=0.1;
    inv_param.maxiter_precondition=10;
    inv_param.verbosity_precondition=get_verbosity_for_quda();
    
    inv_param.omega=1.0;
    
#ifndef DYNAMIC_CLOVER
# warning Please compile quda with DYNAMIC_CLOVER switched on
#endif
    if(exported and csw)
      load_clover_term(&inv_param);
    
    if(multiGrid::checkIfMultiGridAvailableAndRequired(mu))
      {
	inv_param.verbosity=get_verbosity_for_quda();
	
	// coarsening does not support QUDA_MATPC_EVEN_EVEN_ASYMMETRIC
	if(inv_param.matpc_type==QUDA_MATPC_EVEN_EVEN_ASYMMETRIC)
	  inv_param.matpc_type=QUDA_MATPC_EVEN_EVEN;
	
	inv_param.inv_type=QUDA_GCR_INVERTER;
	inv_param.gcrNkrylov=nissa::multiGrid::gcrNkrylov;
	inv_param.inv_type_precondition=QUDA_MG_INVERTER;
	inv_param.schwarz_type=QUDA_ADDITIVE_SCHWARZ;
	inv_param.reliable_delta=nissa::multiGrid::reliable_delta;
	inv_param.reliable_delta_refinement=nissa::multiGrid::reliable_delta_refinement;
	inv_param.precondition_cycle=1;
	inv_param.tol_precondition=1e-1;
	inv_param.maxiter_precondition=1;
	inv_param.gamma_basis=QUDA_CHIRAL_GAMMA_BASIS;
	inv_param.solve_type=QUDA_DIRECT_PC_SOLVE;
	
	inv_param.omega=1.0;
	
	inv_mg_param=inv_param;
	inv_mg_param.preconditioner=nullptr;
	inv_mg_param.inv_type=QUDA_GCR_INVERTER;
	inv_mg_param.inv_type_precondition=QUDA_INVALID_INVERTER;
	inv_mg_param.maxiter=1000;
	inv_mg_param.solve_type=QUDA_DIRECT_SOLVE;
	inv_mg_param.verbosity=get_verbosity_for_quda();
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
	    set_n_vec_batch(quda_mg_param);
	    
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
	    
	    quda_mg_param.verbosity[level]=get_verbosity_for_quda();
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
	    quda_mg_param.coarse_solver_tol[level]=nissa::multiGrid::coarse_solver_tol[level];          //Suggest setting each level to 0.25
	    quda_mg_param.coarse_solver_maxiter[level]=nissa::multiGrid::coarse_solver_maxiter[level];//(level+1==nlevels)?50:100;        //Suggest setting in the range 8-100
	    quda_mg_param.spin_block_size[level]=(level==0)?2:1;  //2 for level 0, and 1 thereafter
	    quda_mg_param.n_vec[level]=(level>=1)?32:24;          //24 or 32 is supported presently
	    quda_mg_param.nu_pre[level]=nissa::multiGrid::nu_pre[level];            //Suggest setting to 0
	    quda_mg_param.nu_post[level]=nissa::multiGrid::nu_post[level];          //Suggest setting to 8
	    
	    //Always set to QUDA_MG_CYCLE_RECURSIVE (this sets the MG cycles to be a K-cycle which is generally superior to a V-cycle for non-Hermitian systems)
	    quda_mg_param.cycle_type[level]=QUDA_MG_CYCLE_RECURSIVE;
	    //Set to QUDA_CUDA_FIELD_LOCATION for all levels
	    quda_mg_param.location[level]=QUDA_CUDA_FIELD_LOCATION;
	    quda_mg_param.setup_location[level]=QUDA_CUDA_FIELD_LOCATION;
	    
	    quda_mg_param.preserve_deflation=QUDA_BOOLEAN_FALSE;
	    quda_mg_param.smoother[level]=(level+1==nlevels)?QUDA_MR_INVERTER:QUDA_CA_GCR_INVERTER;
	    quda_mg_param.smoother_tol[level]=nissa::multiGrid::smoother_tol[level];                 //Suggest setting each level to 0.25
	    quda_mg_param.smoother_schwarz_cycle[level]=1;          //Experimental, set to 1 for each level
	    //Suggest setting to QUDA_DIRECT_PC_SOLVE for all levels
	    quda_mg_param.smoother_solve_type[level]=QUDA_DIRECT_PC_SOLVE;
	    //Experimental, set to QUDA_INVALID_SCHWARZ for each level unless you know what you're doing
	    quda_mg_param.smoother_schwarz_type[level]=QUDA_INVALID_SCHWARZ;
	    quda_mg_param.smoother_halo_precision[level]=QUDA_HALF_PRECISION;
	    
	    // when the Schwarz-alternating smoother is used, this can be set to NO, otherwise it must be YES
	    quda_mg_param.global_reduction[level]=QUDA_BOOLEAN_YES;
	    
	    // set to QUDA_MAT_SOLUTION to inject a full field into coarse grid
	    // set to QUDA_MATPC_SOLUTION to inject single parity field into coarse grid
	    // if we are using an outer even-odd preconditioned solve, then we
	    // use single parity injection into the coarse grid
	    quda_mg_param.coarse_grid_solution_type[level]=
	      QUDA_MATPC_SOLUTION;
	    //(inv_param.solve_type==QUDA_DIRECT_PC_SOLVE?QUDA_MATPC_SOLUTION:QUDA_MAT_SOLUTION);
	    quda_mg_param.omega[level]=nissa::multiGrid::omega[level];  //Set to 0.8-1.0
	    
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
	    if(level+1==nlevels and multiGrid::use_deflated_solver and fabs(inv_param.mu)<multiGrid::max_mass_for_deflation)
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
		mg_eig_param[level].a_min=multiGrid::eig_min;
		mg_eig_param[level].a_max=multiGrid::eig_max;
		
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
    static double storedMu=0;
    static double storedKappa=0;
    static double storedCloverCoeff=0;
    
#ifdef QUDASETUPSTORE
    const char QUDA_DEBUG_EV[]="QUDA_DEBUG_EV";
    const bool doTheStorage=getenv(QUDA_DEBUG_EV)!=nullptr;
    
    SetupID setupId=
      std::make_tuple(export_conf::confTag,export_conf::check_old);
#endif
    
    bool& setup_valid=multiGrid::setup_valid;
    if(not setup_valid)
      {
	master_printf("QUDA multigrid setup not valid\n");
	
#ifdef QUDASETUPSTORE
	const bool canReuseStoredSetup=(qudaSetups.find(setupId)!=qudaSetups.end());
	const auto [tag,fs]=setupId;
	master_printf("CanReuseStoredSetup (%s,%zu,%zu): %s\n",tag.c_str(),fs[0],fs[1],canReuseStoredSetup?"true":"false");
	
	setVerbosity(QUDA_VERBOSE);
	master_printf("VERBOSITY OF QUDA: %d should be %d it might be %d\n",getVerbosity(),QUDA_VERBOSE,QUDA_SUMMARIZE);
	
	if(canReuseStoredSetup)
	  {
	    destroyMultigridQuda(quda_mg_preconditioner);
	    
	    master_printf("mg setup redue:\n");
	    
	    const int& nlevels=multiGrid::nlevels;
	    int b[nlevels];
	    for(int level=0;level<nlevels;level++)
	      {
		int& v=quda_mg_param.num_setup_iter[level];
		b[level]=v;
		v=0;
	      }
	    quda_mg_preconditioner=newMultigridQuda(&quda_mg_param);
	    for(int level=0;level<nlevels;level++)
	      quda_mg_param.num_setup_iter[level]=b[level];
	    qudaSetups[setupId].restore();
	  }
	else
#endif
	  {
	    if(quda_mg_preconditioner!=nullptr)
	      destroyMultigridQuda(quda_mg_preconditioner);
	    
	    quda_mg_preconditioner=newMultigridQuda(&quda_mg_param);
	    
#ifdef READY_TO_DEL
	    createDirac(D,DSloppy,DPre,inv_param,true);
	    
	    using namespace quda;
	    
	    M=(inv_param.inv_type==QUDA_CG_INVERTER or inv_param.inv_type==QUDA_CA_CG_INVERTER)?
	      static_cast<DiracMatrix*>(new DiracMdagM(*D)):
	      static_cast<DiracMatrix*>(new DiracM(*D));
	    MSloppy=(inv_param.inv_type==QUDA_CG_INVERTER or inv_param.inv_type==QUDA_CA_CG_INVERTER)?
	      static_cast<DiracMatrix*>(new DiracMdagM(*DSloppy)):
	      static_cast<DiracMatrix*>(new DiracM(*DSloppy));
	    MPre=(inv_param.inv_type==QUDA_CG_INVERTER or inv_param.inv_type==QUDA_CA_CG_INVERTER)?
	      static_cast<DiracMatrix*>(new DiracMdagM(*DPre)):static_cast<DiracMatrix*>(new DiracM(*DPre));

	    solverParam=new SolverParam(inv_param);
	    
	    const std::string tag="Solver profiler mu="+std::to_string(inv_param.mu);
	    new TimeProfile("ergere");
	    profilers.try_emplace(tag,tag);
	    auto& p=profilers.find(tag)->second;
	    p.TPSTART(QUDA_PROFILE_TOTAL);
	    solver=Solver::create(*solverParam,*M,*MSloppy,*MPre,*MPre,p);
#endif
	    
#ifdef QUDASETUPSTORE
	    if(doTheStorage)
	      qudaSetups[setupId].takeCopy();
#endif
	  }
	
	master_printf("mg setup done!\n");
	
	setup_valid=true;
      }
    else
      if(storedMu!=inv_param.mu or
	 storedKappa!=inv_param.kappa or
	 storedCloverCoeff!=inv_param.clover_coeff)
	{
	  const int& nlevels=
	    multiGrid::nlevels;
	  const double& setup_refresh_tol=
	    multiGrid::setup_refresh_tol;
	  
	  /// Check if tolerance is satisfied, such that only the coarsest level is afjusted
	  const bool tolSatisfied=
	    (fabs(storedMu/(inv_param.mu+1e-300)-1)<setup_refresh_tol) and
	    (fabs(storedKappa/(inv_param.kappa+1e-300)-1)<setup_refresh_tol) and
	    (fabs(storedCloverCoeff/(inv_param.clover_coeff+1e-300)-1)<setup_refresh_tol);
	  master_printf("Tolerance to avoid deep update on csw, mu and kappa change is %lg satisfied: %d\n",setup_refresh_tol,tolSatisfied);
	  
	  int stored_setup_maxiter_refresh[QUDA_MAX_MG_LEVEL];
	  QudaBoolean stored_preserve_deflation;
	  
	  /// Access the quda_mg_param.setup_maxiter_refresh for the asked level
	  auto iR=
	    [](const int& level)->int&
	    {
	      return quda_mg_param.setup_maxiter_refresh[level];
	    };
	  
	  /// Store the parameters, and possibly set it to zero
	  for(int level=0;level<nlevels-1;level++)
	    {
	      stored_setup_maxiter_refresh[level]=iR(level);
	      stored_preserve_deflation=quda_mg_param.preserve_deflation;
	      if(tolSatisfied)
		{
		  iR(level)=0;
		  quda_mg_param.preserve_deflation=QUDA_BOOLEAN_TRUE;
		}
	    }
	  
	  updateMultigridQuda(quda_mg_preconditioner,&quda_mg_param);
	  
	  /// Restore them no matter what
	  for(int level=0;level<nlevels-1;level++)
	    iR(level)=stored_setup_maxiter_refresh[level];
	  quda_mg_param.preserve_deflation=stored_preserve_deflation;
	}
      else
	master_printf("No need to update the multigrid\n");
    
    inv_param.preconditioner=quda_mg_preconditioner;
    
    storedMu=inv_param.mu;
    storedKappa=inv_param.kappa;
    storedCloverCoeff=inv_param.clover_coeff;
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
    for(int ilev=0;ilev<nlev;ilev++)
      printf("transfer_type: %d\n",i.transfer_type[ilev]);
    // printf("use_mma: %d\n",i.use_mma);
    printf("thin_update_only: %d\n",i.thin_update_only);
  }
  
#define PROVIDE_MAYBE_PRINT(X,F)					\
  template <typename T>							\
  void maybe_print_ ## X(const T& i)					\
  {									\
    if constexpr(hasMember_ ## X<T>)					\
      printf(#X": " F "\n",i.X);					\
  }
  
#define PROVIDE_MAYBE_PRINT2(X,F)					\
  template <typename T>							\
  void maybe_print_ ## X(const T& i)					\
  {									\
    if constexpr(hasMember_ ## X<T>)					\
      if constexpr(std::is_array_v<decltype(i.X)>)			\
	printf(#X": " F "\n",i.X[0]);					\
      else								\
	printf(#X": " F "\n",i.X);					\
  }
  
  PROVIDE_MAYBE_PRINT(cl_pad,"%d");
  PROVIDE_MAYBE_PRINT(sp_pad,"%d");
  PROVIDE_MAYBE_PRINT2(true_res,"%lg");
  PROVIDE_MAYBE_PRINT2(true_res_hq,"%lg");
  PROVIDE_MAYBE_PRINT(tune,"%lg");
  
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
    maybe_print_true_res(i);
    maybe_print_true_res_hq(i);
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
    maybe_print_cl_pad(i);
    maybe_print_sp_pad(i);
    printf("iter: %d\n",i.iter);
    printf("gflops: %lg\n",i.gflops);
    printf("secs: %lg\n",i.secs);
    maybe_print_tune(i);
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
  
#undef PROVIDE_MAYBE_PRINT
  
#undef PROVIDE_MAYBE_PRINT2
  
  bool solve_tmD(LxField<spincolor>& sol,
		 const LxField<quad_su3>& conf,
		 const double& kappa,
		 const double& csw,
		 const double& mu,
		 const int& niter,
		 const double& residue,
		 const LxField<spincolor>& source)
  {
    const double export_time=take_time();
    const bool exported=
      export_gauge_conf_to_external_solver(conf);
    
    master_printf("time to export to the conf to quda: %lg s\n",take_time()-export_time);
    
    set_base_inverter_pars();
    
    set_inverter_pars(kappa,csw,mu,niter,residue,exported);
    
#ifdef DEBUG_QUDA
    if(is_master_rank())
      {
	master_printf("--- gauge pars: ---\n");
	printQudaGaugeParam(&gauge_param);
	master_printf("--- inv pars: ---\n");
	printQudaInvertParam(&inv_param);
	
	if(multiGrid::checkIfMultiGridAvailableAndRequired(mu))
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
    
    if(multiGrid::checkIfMultiGridAvailableAndRequired(mu))
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
    crash("reimplement");
    
    // gauge_param.reconstruct=
    //   gauge_param.reconstruct_sloppy=
    //   gauge_param.reconstruct_precondition=
    //   gauge_param.reconstruct_refinement_sloppy=
    //   QUDA_RECONSTRUCT_NO;
    // const double export_time=take_time();
    
    // add_or_rem_stagphases_to_conf(conf);
    // const bool exported=export_gauge_conf_to_external_solver(conf);
    // add_or_rem_stagphases_to_conf(conf);
    // master_printf("time to export (%d) to external library: %lg s\n",exported,take_time()-export_time);
    
    // set_base_inverter_pars();
    
    // inv_param.dslash_type=QUDA_STAGGERED_DSLASH;
    // inv_param.gamma_basis=QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
    // inv_param.matpc_type=QUDA_MATPC_EVEN_EVEN;
    // inv_param.solution_type=QUDA_MAT_SOLUTION;
    
    // inv_param.inv_type=QUDA_CG_INVERTER;
    // inv_param.solve_type=QUDA_DIRECT_PC_SOLVE;
    // inv_param.dagger=QUDA_DAG_NO;
    
    // //minus due to different gamma5 definition
    // inv_param.mass=mass;
    
    // inv_param.verbosity=get_verbosity_for_quda();
    // inv_param.tol=sqrt(residue);
    // inv_param.maxiter=niter;
    // inv_param.Ls=1;
    
    // inv_param.verbosity=get_quda_verbosity();
    
    // remap_nissa_to_quda(color_in,source);
    
    // //invertQuda(color_out,color_in,&inv_param);
    // MatQuda(color_out,color_in,&inv_param);
    
    // master_printf("# QUDA solved in: %i iter / %g secs=%g Gflops\n",inv_param.iter,inv_param.secs,inv_param.gflops/inv_param.secs);
    
    // remap_quda_to_nissa(sol,color_out);
    
    return true;
  }
}
