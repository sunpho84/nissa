#ifndef _FFT_HPP
#define _FFT_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#ifdef USE_CUDA
# include <cufft.h>
#endif

#ifdef USE_FFTW
# include <fftw3.h>
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
	    SpaceTimeLayout STL>
  void fft4d(LxField<T,STL>& f,
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
	constexpr int ncpp=
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
	      
	      /// Number of complexes, needs to be non-const
	      int n=
		glbSize[mu];
	      
	      /// Temporary hack, we need to fuse this step with the previous remap
	      complex *tmp=
		memoryManager<defaultMemorySpace>()->provide<complex>(locVol*ncpp);
	      PAR(0,
		  locd_perp_size_per_dir[mu],
		  CAPTURE(tmp,
			  ncpp,
			  mu,
			  buf),
		  iPerp,
		  {
		    for(int t=0;t<glbSize[mu];t++)
		      for(int icpp=0;icpp<ncpp;icpp++)
			complex_copy(tmp[t+glbSize[mu]*(icpp+ncpp*iPerp)],buf[icpp+ncpp*(t+glbSize[mu]*iPerp)]);
		  });
	      
	      cufftHandle plan;
	      const int stride=1;
	      const int dist=n;
	      decryptFftError(cufftPlanMany(&plan,1,&n,nullptr,stride,dist,nullptr,stride,dist,CUFFT_Z2Z,locd_perp_size_per_dir[mu]*ncpp),"creating the plan");
	      decryptFftError(cufftExecZ2Z(plan,(cufftDoubleComplex*)tmp,(cufftDoubleComplex*)tmp,sign),"executing the transform");
	      DECRYPT_CUDA_ERROR(cudaDeviceSynchronize(),"synchronizing at the end of fft");
	      
	      decryptFftError(cufftDestroy(plan),"destroying the plan");
	      
	      PAR(0,
		  locd_perp_size_per_dir[mu],
		  CAPTURE(tmp,
			  ncpp,
			  mu,
			  buf),
		  iPerp,
		  {
		    for(int t=0;t<glbSize[mu];t++)
		      for(int icpp=0;icpp<ncpp;icpp++)
			complex_copy(buf[icpp+ncpp*(t+glbSize[mu]*iPerp)],tmp[t+glbSize[mu]*(icpp+ncpp*iPerp)]);
		  });
	      
	      memoryManager<defaultMemorySpace>()->release(tmp);
#else
	      
#ifdef USE_FFTW
	      fftw_plan plan=
		fftw_plan_many_dft(1,&glbSize[mu],ncpp,buf,nullptr,ncpp,1,buf,nullptr,ncpp,1,sign,FFTW_ESTIMATE);
	      
	      //makes all the fourier transforms
	      PAR(0,
		  locd_perp_size_per_dir[mu],
		  CAPTURE(plan,buf,mu,ncpp),
		  ioff,
		  {
		    fftw_execute_dft(plan,buf+ioff*glbSize[mu]*ncpp,buf+ioff*glbSize[mu]*ncpp);
		  });
	      
	      fftw_destroy_plan(plan);
#else
	      CRASH("FFTW needed, not found");
#endif
	      
#endif
	      
	      remap_locd_vector_to_lx(fptr,buf,ncpp*sizeof(complex),mu);
	      
	      norm*=glbSize[mu];
	    }
	  
	  if(normalize)
	    f/=norm;
	  
	  memoryManager<defaultMemorySpace>()->release(buf);
      }
  }
}

#endif
