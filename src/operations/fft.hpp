#ifndef _FFT_HPP
#define _FFT_HPP

#include <fftw3.h>

#ifdef USE_CUDA
# include <cufft.h>
#endif

#include "base/field.hpp"
#include "operations/remap_vector.hpp"

namespace nissa
{
  enum{FFT_NO_NORMALIZE=0,FFT_NORMALIZE=1};
#define FFT_PLUS +1.0
#define FFT_MINUS -1.0
  
  /// Fourier transform of the field
  template <typename T,
	    SpaceTimeLayout STL,
	    MemorySpace MS>
  void fft4d(LxField<T,STL,MS>& f,
	     const double& sign,
	     const bool& normalize)
  {
    static_assert(LxField<T>::nInternalDegs%2==0,"Needs to be able to interpret the field as complex");
    
    if constexpr(STL!=SpaceTimeLayout::CPU)
      {
	LxField<T,SpaceTimeLayout::CPU> tmp("tmp");
	tmp=f;
	MASTER_PRINTF("Making a temporary copy of the field changing layout to be 4D-fourier-transformed\n");
	
	fft4d(tmp,sign,normalize);
	
	f=tmp;
      }
    else
      {
	const int ncpp=
	  LxField<T>::nInternalDegs/2;
	
	double norm=1;
	
	//allocate buffer
	complex *buf=
	  memoryManager<defaultMemorySpace>()->provide<complex>(max_locd_size*ncpp);
	
	//transpose each dir in turn and take fft
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      VERBOSITY_LV2_MASTER_PRINTF("FFT-ing dimension %d\n",mu);
	      
	      auto* fptr=
		f.template getPtr<defaultMemorySpace>();
	      
	      remap_lx_vector_to_locd(buf,fptr,ncpp*sizeof(complex),mu);
	      
#ifdef USE_CUDA
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
	      
	      /// Number of complexes, needs to be non-const
	      int n=
		glbSize[mu];
	      
	      /// Extension of the fft
	      const int64_t nFft=
		locd_perp_size_per_dir[mu]/ncpp;
	      
	      decryptFftError(cufftPlanMany(&plan,1,&n,&n,1,n,&n,1,n,CUFFT_Z2Z,nFft),"creating the plan");
	      decryptFftError(cufftExecZ2Z(plan,b,b,sign),"executing the transform");
	      DECRYPT_CUDA_ERROR(cudaDeviceSynchronize(),"synchronizing at the end of fft");
	      
	      decryptFftError(cufftDestroy(plan),"destroying the plan");
#else
	      fftw_plan plan=
		fftw_plan_many_dft(1,&glbSize[mu],ncpp,buf,nullptr,ncpp,1,buf,nullptr,ncpp,1,sign,FFTW_ESTIMATE);
	      
	      //makes all the fourier transform
	      PAR(0,
		  locd_perp_size_per_dir[mu],
		  CAPTURE(plan,buf,mu,ncpp),
		  ioff,
		  {
		    fftw_execute_dft(plan,buf+ioff*glbSize[mu]*ncpp,buf+ioff*glbSize[mu]*ncpp);
		  });
	      
	      fftw_destroy_plan(plan);
#endif
	      
	      remap_locd_vector_to_lx(fptr,buf,ncpp*sizeof(complex),mu);
	      
	      norm*=glbSize[mu];
	    }
	  
	  f/=norm;
	  
	  memoryManager<defaultMemorySpace>()->release(buf);
      }
  }
}

#endif
