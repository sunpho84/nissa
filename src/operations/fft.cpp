#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <fftw3.h>

#ifdef USE_CUDA
# include <cufft.h>
#endif

#define EXTERN_FFT
# include <operations/fft.hpp>

#include <routines/ios.hpp>

namespace nissa
{
  void initFftw()
  {
    if(fftwInitialized)
      crash("cannot initialize fftw twice");
    
    master_printf("Initializing fftw\n");
    
#ifdef USE_FFTW_THREADED
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif
    
    fftwInitialized=true;
  }
  
  void fftwFinalize()
  {
    if(not fftwInitialized)
      crash("fftw not initialized");
    
    master_printf("Finalizing fftw\n");
    
#ifdef USE_FFTW_THREADED
    fftw_cleanup_threads();
#endif
    
    fftwInitialized=false;
  }
  
  void fftExecUsingFftw(void* buf,
			const int& n,
			const int& sign,
			const int& nFft)
  {
    const double startTime=
      take_time();
    
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
    
    master_printf("fft executed on CPU "
#ifdef USE_FFTW_THREADED
		  "with fftw-threads "
#endif
		  "in %lg s\n",take_time()-startTime);
  }
  
#ifdef USE_CUDA
  void fftExecUsingCuFFT(void* buf,
			 int n,
			 const int& sign,
			 const int& nFft)
  {
    const double startTime=
      take_time();

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
	    crash("fft crashed at stage %s with error %s",stage,name);
	
	crash("fft crashed at stage %s with unknown error %d %d",stage,res,CUFFT_SUCCESS);
      };
    
    cufftHandle plan;
    decryptFftError(cufftPlanMany(&plan,1,&n,&n,1,n,&n,1,n,CUFFT_Z2Z,nFft),"creating the plan");
    decryptFftError(cufftExecZ2Z(plan,b,b,sign),"executing the transform");
    decrypt_cuda_error(cudaDeviceSynchronize(),"synchronizing at the end of fft");
    
    decryptFftError(cufftDestroy(plan),"destroying the plan");
    
    master_printf("fft executed on GPU in %lg s\n",take_time()-startTime);
  }
#endif
}
