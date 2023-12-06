#ifndef _CUDA_TEST_HPP
#define _CUDA_TEST_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <geometry/geometry_eo.hpp>
#include <new_types/su3.hpp>
#include <routines/ios.hpp>

#ifdef USE_CUDA

namespace nissa
{
  /// Complex number with base type T
  template <typename T>
  struct Compl
  {
    /// Real part
    T re;
    
    /// Imaginary part
    T im;
    
    /// Set from a complex
    __host__ __device__ void set(const complex &o)
    {
      re=o[RE];
      im=o[IM];
    }
    
    /// Assign
    template <typename O>
    __host__ __device__ Compl& operator=(const O& oth)
    {
      re=(O)oth;
      im=T{0};
      
      return *this;
    }
    
    /// Take the conjugate
    __host__ __device__ Compl conj() const
    {
      return {re,-im};
    }
    
    /// Product
    template <typename O>
    __host__ __device__ Compl<decltype(O{}*T{})> operator*(const Compl<O>& oth) const
    {
      return {re*oth.re-im*oth.im,re*oth.im+im*oth.re};
    }
    
    /// Summassign
    template <typename O>
    __host__ __device__ Compl<decltype(O{}+T{})> operator+=(const Compl<O>& oth)
    {
      re+=oth.re;
      im+=oth.im;
      
      return *this;
    }
  };
  
  namespace gpu
  {
    /// Allocate on gpu
    template <typename T>
    void gpu_alloc(T*& data,const int64_t n)
    {
      const int64_t tot=n*sizeof(T);
      
      decrypt_cuda_error(cudaMalloc(&data,tot),"Allocating %ld bytes",tot);
    }
    
    /// Free on gpu
    template <typename T>
    void gpu_free(T*& data)
    {
      decrypt_cuda_error(cudaFree(data),"Freeing");
      data=NULL;
    }
    
    /// Move to gpu
    template <typename T>
    void cpu_to_gpu(T* out,const T* in,const int64_t n)
    {
      const int64_t tot=n*sizeof(T);
      decrypt_cuda_error(cudaMemcpy(out,in,tot,cudaMemcpyHostToDevice),"Copying %ld bytes from cpu to gpu",tot);
    }
    
    /// Move to cpu
    template <typename T>
    void gpu_to_cpu(T* out,const T* in,const int64_t n)
    {
      const int64_t tot=n*sizeof(T);
      decrypt_cuda_error(cudaMemcpy(out,in,n*sizeof(T),cudaMemcpyDeviceToHost),"Copying %ld bytes from cpu to gpu",tot);
    }
    
    /// Color on a gpu
    template <typename T>
    class gpu_color
    {
      Compl<T>* data;
      
    public:
      
      /// Index of the data
      __host__ __device__ int64_t idx(const int icol,const int64_t ivol_eo) const
      {
	return ivol_eo+locVolh*icol;
	//return icol+NCOL*ivol_eo;
      }
      
      /// Constant access to the data
      __device__ const Compl<T>& operator()(const int icol,const int64_t ivol_eo) const
      {
	return data[idx(icol,ivol_eo)];
      }
      
      /// Access to the data
      __device__ Compl<T>& operator()(const int icol,const int64_t ivol_eo)
      {
	return data[idx(icol,ivol_eo)];
      }
      
      /// Size
      const int64_t n;
      
      /// Default constructor
      gpu_color() : n(NCOL*(locVolh+bord_volh))
      {
      }
      
      /// Allocate
      void alloc()
      {
	gpu_alloc(data,n);
      }
      
      /// Deallocate
      void dealloc()
      {
	gpu_free(data);
      }
      
      /// Import on the gpu the passed cpu vector
      void import_on_gpu(const color *in)
      {
	std::unique_ptr<Compl<T>[]> buf(new Compl<T>[n]);
	
	for(int64_t ivol_eo=0;ivol_eo<locVolh;ivol_eo++)
	  for(int icol=0;icol<NCOL;icol++)
	    {
	      const int64_t iout=idx(icol,ivol_eo);
	      const complex& temp=in[ivol_eo][icol];
	      buf[iout].set(temp);
	    }
	
	master_printf("Preparing transfer\n");
	
	cpu_to_gpu(data,&buf[0],n);
	
	master_printf("Transferred color\n");
      }
      
      /// Export to the passed cpu vector
      void export_to_cpu(color *out) const
      {
	std::unique_ptr<Compl<T>[]> buf(new Compl<T>[n]);
	
	gpu_to_cpu(&buf[0],data,n);
	
	for(int icol=0;icol<NCOL;icol++)
	  for(int64_t ivol_eo=0;ivol_eo<locVolh;ivol_eo++)
	    {
	      const int64_t iin=idx(icol,ivol_eo);
	      complex& temp=out[ivol_eo][icol];
	      temp[RE]=buf[iin].re;
	      temp[IM]=buf[iin].im;
	    }
      }
    };
    
    /////////////////////////////////////////////////////////////////
    
    /// Conf on the gpu
    template <typename T>
    class gpu_links
    {
      /// Data
      Compl<T>* data;
      
    public:
      
      /// Index
      __host__ __device__ int64_t idx(const int mu,const int icol1,const int icol2,const int ieo,const int64_t ivol_eo) const
      {
	return ivol_eo+locVolh*(ieo+2*(icol2+NCOL*(icol1+NCOL*mu)));
	//return ivol_eo+locVolh*(ieo+2*(mu+NDIM*(icol2+NCOL*icol1)));
	//return icol2+NCOL*(icol1+NCOL*(mu+NDIM*(ivol_eo+locVolh*ieo)));
      }
      
      /// Access to the data
      __device__ Compl<T>& operator()(const int mu,const int icol1,const int icol2,const int ieo,const int64_t ivol_eo)
      {
	return data[idx(mu,icol1,icol2,ieo,ivol_eo)];
      }
      
      /// Constant access to the data
      __device__ const Compl<T>& operator()(const int mu,const int icol1,const int icol2,const int ieo,const int64_t ivol_eo) const
      {
	return data[idx(mu,icol1,icol2,ieo,ivol_eo)];
      }
      
      /// Size
      const int64_t n;
      
      /// Default constructor
      gpu_links() : n(NDIM*NCOL*NCOL*2*(locVolh+bord_volh))
      {
      }
      
      /// Allocate
      void alloc()
      {
	gpu_alloc(data,n);
      }
      
      /// Deallocate
      void dealloc()
      {
	gpu_free(data);
      }
      
      /// Copy on the gpu the passed conf
      void import_on_gpu(const eo_ptr<quad_su3> in)
      {
	std::unique_ptr<Compl<T>[]> buf(new Compl<T>[n]);
	
	for(int ieo=0;ieo<2;ieo++)
	  for(int64_t ivol_eo=0;ivol_eo<locVolh;ivol_eo++)
	    for(int mu=0;mu<NDIM;mu++)
	      for(int icol1=0;icol1<NCOL;icol1++)
		for(int icol2=0;icol2<NCOL;icol2++)
		  {
		    const int64_t iout=idx(mu,icol1,icol2,ieo,ivol_eo);
		    const complex& temp=in[ieo][ivol_eo][mu][icol1][icol2];
		    buf[iout].set(temp);
		  }
	
	master_printf("Preparing transfer\n");
	
	cpu_to_gpu(data,&buf[0],n);
	
	master_printf("Transferred conf\n");
      }
      
      /// Copy to the passed cpu vector
      void export_to_cpu(eo_ptr<quad_su3> out) const
      {
	std::unique_ptr<Compl<T>[]> buf(new Compl<T>[n]);
	
	cpu_to_gpu(&buf[0],data,n);
	
	for(int mu=0;mu<NDIM;mu++)
	  for(int icol1=0;icol1<NCOL;icol1++)
	    for(int icol2=0;icol2<NCOL;icol2++)
	      for(int ieo=0;ieo<2;ieo++)
		for(int64_t ivol_eo=0;ivol_eo<locVolh;ivol_eo++)
		  {
		    const int64_t iin=idx(mu,icol1,icol2,ieo,ivol_eo);
		    complex& temp=out[ieo][ivol_eo][mu][icol1][icol2];
		    temp[RE]=buf[iin].real();
		    temp[IM]=buf[iin].imag();
		  }
      }
    };
    
    template <int PAR,
	      typename T>
    __global__ void Doe_or_Deo(gpu_color<T> out,const gpu_links<T> conf,const gpu_color<T> in)
    {
      
      const int64_t ivol_out=blockIdx.x*blockDim.x+threadIdx.x;
      if(ivol_out<locVolh)
	{
	  for(int ic=0;ic<NCOL;ic++)
	    out(ic,ivol_out)=0.0;
	  
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      const int64_t ivol_up_in=
		loceo_neighup[PAR][ivol_out][mu];
	      
	      for(int ic1=0;ic1<NCOL;ic1++)
	  	for(int ic2=0;ic2<NCOL;ic2++)
	  	  out(ic1,ivol_out)+=conf(mu,ic1,ic2,PAR,ivol_out)*in(ic2,ivol_up_in);
		  // out(ic1,ivol_out).summ_the_prod(conf(mu,ic1,ic2,PAR,ivol_out),in(ic2,ivol_up_in));
	      
	      const int64_t ivol_dw_in=
		loceo_neighdw[PAR][ivol_out][mu];
	      
	      for(int ic1=0;ic1<NCOL;ic1++)
	  	for(int ic2=0;ic2<NCOL;ic2++)
		  out(ic1,ivol_out)+=conf(mu,ic1,ic2,!PAR,ivol_dw_in).conj()*in(ic2,ivol_dw_in);
	      //out(ic1,ivol_out).summ_the_prod(conf(mu,ic1,ic2,!PAR,ivol_dw_in),in(ic2,ivol_dw_in));
	    }
	}
      // else
      // 	printf("%d not running\n",ivol_out);
    }
    
    template <typename T>
    void cuda_test(color *_out,eo_ptr<quad_su3> _conf,color *_in)
    {
      gpu_color<T> in;
      in.alloc();
      in.import_on_gpu(_in);
      
      gpu_color<T> temp;
      temp.alloc();
      
      gpu_links<T> conf;
      conf.alloc();
      conf.import_on_gpu(_conf);
      
      gpu_color<T> out;
      out.alloc();
      
      for(int ngpu_threads=2;ngpu_threads<1024;ngpu_threads*=2)
	{
	  master_printf("nthreads: %d\n",ngpu_threads);
      const dim3 block_dimension(ngpu_threads);
      const dim3 grid_dimension((locVolh+ngpu_threads)/ngpu_threads);
      
      double init=take_time();
      int n=100;
      
      for(int i=0;i<n;i++)
	{
	  Doe_or_Deo<EVN><<<grid_dimension,block_dimension>>>(temp,conf,in);
	  cudaDeviceSynchronize();
	  Doe_or_Deo<ODD><<<grid_dimension,block_dimension>>>(out,conf,temp);
	  cudaDeviceSynchronize();
	}
      
      double end=take_time();
      
      double each=(end-init)/n;
      const int nflops_per_site=8*8*9*2;
      master_printf("Time for the improved operator: %lg s, per site: %lg s, Gflops: %lg\n",each,each/locVolh,nflops_per_site/each*locVolh*1e-9);
	}
      out.export_to_cpu(_out);
      
      out.dealloc();
      conf.dealloc();
      temp.dealloc();
      in.dealloc();
    }
  }
}

#else

namespace nissa
{
  template <typename T>
  void cuda_test(color *_out,eo_ptr<quad_su3> _conf,color *_in)
  {
  }
}

#endif

#endif
