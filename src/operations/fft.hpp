#ifndef _FFT_HPP
#define _FFT_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#ifndef EXTERN_FFT
# define EXTERN_FFT extern
# define INITIALIZE_FFT_TO(ARGS...)
#else
# define INITIALIZE_FFT_TO(ARGS...) ARGS
#endif

#include <expr/conj.hpp>
#include <expr/field.hpp>

namespace nissa
{
  EXTERN_FFT bool fftwInitialized INITIALIZE_FFT_TO({false});
  
  void initFftw();
  
  void fftwFinalize();
  
  void fftExecUsingFftw(void* buf,const int& n,const int& sign,const int& nFft);
  
  void fftExecUsingCuFFT(void* buf,int n,const int& sign,const int& nFft);
  
  /// Performs fft of in  with given sign putting the result in out
  template <DerivedFromNode Out,
	    DerivedFromNode In>
  void fft(Out&& out,
	   const int& sign,
	   const In& in)
  {
    /// Components different from those related to fft
    using OthComps=
      TupleFilterAllTypes<typename In::Comps,CompsList<LocLxSite,ComplId>>;
    
    auto f=
      [sign]<DerivedFromNode B,
	     Dir D>(B& buf,
		    const std::integral_constant<Dir,D>&)
      {
	const int nCompl=glbSizes[D()];
	const int nFft=buf.nElements/nCompl/2;
	
#ifdef USE_CUDA
	if constexpr(B::execSpace==execOnGPU)
	  fftExecUsingCuFFT(buf.storage,nCompl,sign,nFft);
	else
#endif
	  fftExecUsingFftw(buf.storage,nCompl,sign,nFft);
	
	master_printf("%s\n",demangle(typeid(B).name()).c_str());
	
	master_printf("FFTing on Dir %d nFft=%d\n",D(),nFft);
      };
    
    cycleOnAllLocalDirections<OthComps,CompsList<ComplId>>(std::forward<Out>(out),in,f,in.getDynamicSizes());
  }
  
  template <DerivedFromNode In>
  auto fft(const int& sign,
	   const In& in)
  {
    In out;
    
    fft(out,sign,in);
    
    return out;
  }
}

#undef EXTERN_FFT
#undef INITIALIZE_FFT_TO

#endif
