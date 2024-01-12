#ifndef _FFT_HPP
#define _FFT_HPP

#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <fftw3.h>

#ifdef ENABLE_DEVICE_CODE
# include <cufft.h>
#endif

#ifdef USE_OPENMP
# include <omp.h>
#endif

#include <expr/conj.hpp>
#include <expr/field.hpp>

namespace nissa
{
  inline bool fftwInitialized{false};
  
  void fftExecUsingFftw(void* buf,const int& n,const int& sign,const int& nFft);
  
  void fftExecUsingCuFFT(void* buf,int n,const int& sign,const int& nFft);
  
  /// Performs fft of in  with given sign putting the result in out
  template <DerivedFromNode Out,
	    DerivedFromNode In>
  void fft(Out&& out,
	   const int& sign,
	   const In& in)
  {
    if constexpr(tupleHasType<typename In::Comps,ComplId>)
      {
	/// Components different from those related to fft
	using OthComps=
	  TupleFilterAllTypes<typename In::Comps,CompsList<LocLxSite,ComplId>>;
	
	/// Function to exec at each loop iteration
	auto f=
	  [sign]<DerivedFromNode B,
		 Dir D>(B& buf,
			const std::integral_constant<Dir,D>&)
	  {
	    /// Number of complexes
	    const int64_t nCompl=lat->getGlbSizes()[D]();
	    
	    /// Extension of the fft
	    const int64_t nFft=buf.nElements/nCompl/2;
	    
#ifdef ENABLE_DEVICE_CODE
	    if constexpr(B::execSpace==execOnGPU)
	      fftExecUsingCuFFT(buf.storage,nCompl,sign,nFft);
	    else
#endif
	      fftExecUsingFftw(buf.storage,nCompl,sign,nFft);
	    
	    masterPrintf("%s\n",demangle(typeid(B).name()).c_str());
	    
	    masterPrintf("FFTing on Dir %d nFft=%ld\n",D(),nFft);
	  };
	
	cycleOnAllLocalDirections<OthComps,CompsList<ComplId>>(std::forward<Out>(out),in,f,in.getDynamicSizes());
      }
    else
      fft(std::forward<Out>(out),sign,in*complOne<>);
  }
  
  /// Takes fft and output the result to a new object of type Out
  template <DerivedFromNode Out,
	    DerivedFromNode In>
  Out fft(const int& sign,
	  const In& in)
  {
    /// Result
    Out out;
    
    fft(out,sign,in);
    
    return out;
  }
  
  /// Initializes the fftw library
  inline void initFftw()
  {
    if(fftwInitialized)
      CRASH("cannot initialize fftw twice");
    
    masterPrintf("Initializing fftw\n");
    
#ifdef USE_FFTW_THREADED
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif
    
    fftwInitialized=true;
  }
  
  /// Finalizes the fftw library
  inline void fftwFinalize()
  {
    if(not fftwInitialized)
      CRASH("fftw not initialized");
    
    masterPrintf("Finalizing fftw\n");
    
#ifdef USE_FFTW_THREADED
    fftw_cleanup_threads();
#endif
    
    fftwInitialized=false;
  }
  
  /// Exec the fft using libfftw
  inline void fftExecUsingFftw(void* buf,
			       const int& n,
			       const int& sign,
			       const int& nFft)
  {
    const double startTime=
      takeTime();
    
    auto b=
      (fftw_complex*)buf;
    
#ifdef USE_FFTW_THREADED
    fftw_plan plan=
      fftw_plan_many_dft(1,&n,nFft,b,&n,1,n,b,&n,1,n,sign,FFTW_ESTIMATE);
    
    fftw_execute_dft(plan,b,b);
    
    fftw_destroy_plan(plan);
#else
    fftw_plan plan=
      fftw_plan_many_dft(1,&n,1,b,nullptr,1,1,b,nullptr,1,1,sign,FFTW_ESTIMATE);
    
    HOST_PARALLEL_LOOP(0,
		       nFft,
		       CAPTURE(b,
			       plan,
			       n),
		       ioff,
		       {
			 fftw_execute_dft(plan,b+ioff*n,b+ioff*n);
		       });
    
    fftw_destroy_plan(plan);
#endif
    
    masterPrintf("fft executed on CPU "
#ifdef USE_FFTW_THREADED
		  "with fftw-threads "
#endif
		  "in %lg s\n",takeTime()-startTime);
  }
  
#ifdef ENABLE_DEVICE_CODE
  /// Performs the fft using cuda
  inline void fftExecUsingCuFFT(void* buf,
				int n,
				const int& sign,
				const int& nFft)
  {
    const double startTime=
      takeTime();
    
    auto b=
      (cufftDoubleComplex*)buf;
    
    auto decryptFftError=
      [](const cufftResult& res,
	 const char* stage)
      {
	if(res==CUFFT_SUCCESS)
	  return ;
	
	for(const auto& [err,name] :
#define E(A)					\
	      std::make_tuple(CUFFT_ ## A,#A)
	      {E(INVALID_PLAN),E(ALLOC_FAILED),E(INVALID_VALUE),E(INTERNAL_ERROR),E(SETUP_FAILED),E(INVALID_SIZE)})
#undef E
	  if(err==res)
	    CRASH("fft crashed at stage %s with error %s",stage,name);
	
	CRASH("fft crashed at stage %s with unknown error %d %d",stage,res,CUFFT_SUCCESS);
      };
    
    cufftHandle plan;
    decryptFftError(cufftPlanMany(&plan,1,&n,&n,1,n,&n,1,n,CUFFT_Z2Z,nFft),"creating the plan");
    decryptFftError(cufftExecZ2Z(plan,b,b,sign),"executing the transform");
    decryptCudaError(cudaDeviceSynchronize(),"synchronizing at the end of fft");
    
    decryptFftError(cufftDestroy(plan),"destroying the plan");
    
    masterPrintf("fft executed on GPU in %lg s\n",takeTime()-startTime);
  }
#endif
}

#endif
