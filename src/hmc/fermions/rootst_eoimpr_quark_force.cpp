#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <complex>
#include <memory>

#include "base/bench.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "hmc/backfield.hpp"
#include "inverters/staggered/cgm_invert_stD2ee_m2.hpp"
#include "new_types/su3.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  template <typename T>
  struct Compl
  {
    T re;
    T im;
    
    __device__ Compl& operator=(const T& oth)
    {
      re=oth;
      im=T{0};
      
      return *this;
    }
    
    template <typename O>
    __device__ Compl<decltype(O{}*T{})> operator*(const Compl<O>& oth) const
    {
      return {re*oth.re-im*oth.im,re*oth.im+im*oth.re};
    }
    
    template <typename O>
    __device__ Compl<decltype(O{}+T{})> operator+=(const Compl<O>& oth)
    {
      re+=oth.re;
      im+=oth.im;
      
      return *this;
    }
  };
  
  namespace gpu
  {
    template <typename T>
    void gpu_alloc(T*& data,const int64_t n)
    {
      cudaMalloc(&data,n*sizeof(T));
    }
    
    template <typename T>
    void gpu_free(T& data)
    {
      cudaFree(&data);
    }
    
    template <typename T>
    void cpu_to_gpu(T* out,const T* in,const int64_t n)
    {
      cudaMemcpy(out,in,n*sizeof(T),cudaMemcpyHostToDevice);
    }
    
    template <typename T>
    void gpu_to_cpu(T* out,const T* in,const int64_t n)
    {
      cudaMemcpy(out,in,n*sizeof(T),cudaMemcpyDeviceToHost);
    }
    
    template <typename T>
    class gpu_color
    {
      Compl<T>* data;
      
    public:
      
      __host__ __device__ int64_t idx(const int icol,const int64_t ivol_eo) const
      {
	return ivol_eo+loc_volh*icol;
      }
      
      __device__ const Compl<T>& operator()(const int icol,const int64_t ivol_eo) const
      {
	return data[idx(icol,ivol_eo)];
      }
      
      __device__ Compl<T>& operator()(const int icol,const int64_t ivol_eo)
      {
	return const_cast<Compl<T>&>(const_cast<gpu_color&>(*this)(icol,ivol_eo));
      }
      
      const int64_t n;
      
      gpu_color(const color *in=nullptr): n(NCOL*loc_volh)
      {
	gpu_alloc(data,n);
	if(in!=nullptr)
	  import_on_gpu(in);
      }
      
      ~gpu_color()
      {
	gpu_free(data);
      }
      
      void import_on_gpu(const color *in)
      {
	std::unique_ptr<Compl<T>[]> buf(new Compl<T>[n]);
	
	for(int64_t ivol_eo=0;ivol_eo<loc_volh;ivol_eo++)
	  for(int icol=0;icol<NCOL;icol++)
	    {
	      int64_t iout=idx(icol,ivol_eo);
	      const complex& temp=in[ivol_eo][icol];
	      buf[iout]={temp[RE],temp[IM]};
	    }
	
	cpu_to_gpu(data,&buf[0],n);
      }
      
      void export_to_cpu(color *out) const
      {
	std::unique_ptr<Compl<T>[]> buf(new Compl<T>[n]);
	
	cpu_to_gpu(&buf[0],data,n);
	
	for(int icol=0;icol<NCOL;icol++)
	  for(int64_t ivol_eo=0;ivol_eo<loc_volh;ivol_eo++)
	    {
	      const int64_t iin=idx(icol,ivol_eo);
	      complex& temp=out[ivol_eo][icol];
	      temp[RE]=buf[iin].re;
	      temp[IM]=buf[iin].im;
	    }
      }
    };
    
    /////////////////////////////////////////////////////////////////
    
    template <typename T>
    class gpu_links
    {
      Compl<T>* data;
      
    public:
      
      __host__ __device__ int64_t idx(const int mu,const int icol1,const int icol2,const int ieo,const int64_t ivol_eo) const
      {
	return ivol_eo+loc_volh*(ieo+2*(icol2+NCOL*(icol1+NCOL*mu)));
      }
      
      __device__ Compl<T>& operator()(const int mu,const int icol1,const int icol2,const int ieo,const int64_t ivol_eo)
      {
	const_cast<Compl<T>&>(const_cast<gpu_links&>(*this)(mu,icol1,icol2,ieo,ivol_eo));
      }
      
      __device__ const Compl<T>& operator()(const int mu,const int icol1,const int icol2,const int ieo,const int64_t ivol_eo) const
      {
	return data[idx(mu,icol1,icol2,ieo,ivol_eo)];
      }
      
      const int64_t n;
      
      gpu_links(const eo_ptr<quad_su3> in) : n(NDIM*NCOL*NCOL*2*loc_volh)
      {
	gpu_alloc(data,n);
	if(in!=nullptr)
	  import_on_gpu(in);
      }
      
      ~gpu_links()
      {
	gpu_free(data);
      }
      
      void import_on_gpu(const eo_ptr<quad_su3> in)
      {
	std::unique_ptr<Compl<T>[]> buf(new Compl<T>[n]);
	
	for(int ieo=0;ieo<2;ieo++)
	  for(int64_t ivol_eo=0;ivol_eo<loc_volh;ivol_eo++)
	    for(int mu=0;mu<NDIM;mu++)
	      for(int icol1=0;icol1<NCOL;icol1++)
		for(int icol2=0;icol2<NCOL;icol2++)
		  {
		    const int64_t iout=idx(mu,icol1,icol2,ieo,ivol_eo);
		    const complex& temp=in[ieo][ivol_eo][mu][icol1][icol2];
		    buf[iout]={temp[RE],temp[IM]};
		  }
	
	cpu_to_gpu(data,&buf[0],n);
      }
      
      void export_to_cpu(eo_ptr<quad_su3> out) const
      {
	std::unique_ptr<Compl<T>[]> buf(new Compl<T>[n]);
	
	cpu_to_gpu(&buf[0],data,n);
	
	for(int mu=0;mu<NDIM;mu++)
	  for(int icol1=0;icol1<NCOL;icol1++)
	    for(int icol2=0;icol2<NCOL;icol2++)
	      for(int ieo=0;ieo<2;ieo++)
		for(int64_t ivol_eo=0;ivol_eo<loc_volh;ivol_eo++)
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
    __global__ void Doe_or_Deo(gpu_color<T>& out,const gpu_links<T>& conf,const gpu_color<T>& in)
    {
      const int64_t ivol_out=blockIdx.x*blockDim.x+threadIdx.x;
      if(ivol_out<loc_volh)
	{
	  for(int ic=0;ic<NCOL;ic++)
	    out(ic,ivol_out)=0.0;
	  
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      const int64_t ivol_up_in=loceo_neighup[PAR][ivol_out][mu];
	      
	      for(int ic1=0;ic1<NCOL;ic1++)
		for(int ic2=0;ic2<NCOL;ic2++)
		  out(ic1,ivol_out)+=conf(mu,ic1,ic2,PAR,ivol_out)*in(ic2,ivol_up_in);
	      
	      const int64_t ivol_dw_in=loceo_neighdw[PAR][ivol_out][mu];
	      
	      for(int ic1=0;ic1<NCOL;ic1++)
		for(int ic2=0;ic2<NCOL;ic2++)
		  out(ic1,ivol_out)+=conf(mu,ic1,ic2,!PAR,ivol_dw_in)*in(ic2,ivol_dw_in);
	    }
	}
    }
    
    template <typename T>
    void operator_test(color *_out,eo_ptr<quad_su3> _conf,color *_in)
    {
      gpu_color<T> in(_in);
      
      gpu_color<T> temp;
      
      gpu_links<T> conf(_conf);
      
      gpu_color<T> out;
      
      const int ngpu_threads=128;
      const dim3 block_dimension(ngpu_threads);
      const dim3 grid_dimension((loc_volh+ngpu_threads)/ngpu_threads);
      
      double init=take_time();
      int n=100;
      
      if(0)
	for(int i=0;i<n;i++)
	{
	  Doe_or_Deo<EVN><<<grid_dimension,block_dimension>>>(temp,conf,in);
	  Doe_or_Deo<ODD><<<grid_dimension,block_dimension>>>(out,conf,temp);
	}
      
      double end=take_time();
      
      double each=(end-init)/n;
      master_printf("Time for the improved operator: %lg s\n",each);
      
      out.export_to_cpu(_out);
    }
  }
  
  //Compute the fermionic force the rooted staggered eoprec improved theory.
  //Of the result still need to be taken the TA and product with U
  //The approximation need to be already scaled, and must contain physical mass term
  THREADABLE_FUNCTION_8ARG(summ_the_rootst_eoimpr_quark_force, eo_ptr<quad_su3>,F, double,charge, eo_ptr<quad_su3>,eo_conf, color*,pf, int,quantization, eo_ptr<quad_u1>,u1b, rat_approx_t*,appr, double,residue)
  {
    GET_THREAD_ID();
    
    const int nterms=appr->degree();
    
    START_TIMING(quark_force_over_time,nquark_force_over);
    
    //allocate each terms of the expansion
    color **v_o=nissa_malloc("v_o",nterms,color*),**chi_e=nissa_malloc("chi_e",nterms,color*);
    for(int iterm=0;iterm<nterms;iterm++)
      {
	v_o[iterm]=nissa_malloc("v_o",loc_volh+bord_volh,color);
	chi_e[iterm]=nissa_malloc("chi_e",loc_volh+bord_volh,color);
      }
    
    gpu::operator_test<double>(v_o[0],eo_conf,pf);
    
    //add the background fields
    add_backfield_with_stagphases_to_conf(eo_conf,u1b);
    
    //invert the various terms
    STOP_TIMING(quark_force_over_time);
    inv_stD2ee_m2_cgm_run_hm_up_to_comm_prec(chi_e,eo_conf,appr->poles.data(),nterms,1000000,residue,pf);
    UNPAUSE_TIMING(quark_force_over_time);
    
    ////////////////////
    
    //summ all the terms performing appropriate elaboration
    //possible improvement by communicating more borders together
    for(int iterm=0;iterm<nterms;iterm++) apply_stDoe(v_o[iterm],eo_conf,chi_e[iterm]);
    
    //remove the background fields
    rem_backfield_with_stagphases_from_conf(eo_conf,u1b);
    
    //communicate borders of v_o (could be improved...)
    for(int iterm=0;iterm<nterms;iterm++) communicate_od_color_borders(v_o[iterm]);
    
    //conclude the calculation of the fermionic force
    for(int iterm=0;iterm<nterms;iterm++)
      {
	const double weight=appr->weights[iterm];
	
	NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	for(int mu=0;mu<NDIM;mu++)
	  for(int ic1=0;ic1<NCOL;ic1++)
	    for(int ic2=0;ic2<NCOL;ic2++)
	      {
		complex temp1,temp2;
		
		//this is for ieo=EVN
		unsafe_complex_conj2_prod(temp1,v_o[iterm][loceo_neighup[EVN][ieo][mu]][ic1],chi_e[iterm][ieo][ic2]);
		unsafe_complex_prod(temp2,temp1,u1b[EVN][ieo][mu]);
		complex_summ_the_prod_double(F[EVN][ieo][mu][ic1][ic2],temp2,weight*get_stagphase_of_lx(loclx_of_loceo[EVN][ieo],mu));
		
		//this is for ieo=ODD
		unsafe_complex_conj2_prod(temp1,chi_e[iterm][loceo_neighup[ODD][ieo][mu]][ic1],v_o[iterm][ieo][ic2]);
		unsafe_complex_prod(temp2,temp1,u1b[ODD][ieo][mu]);
		complex_subt_the_prod_double(F[ODD][ieo][mu][ic1][ic2],temp2,weight*get_stagphase_of_lx(loclx_of_loceo[ODD][ieo],mu));
	      }
	NISSA_PARALLEL_LOOP_END;
      }
    
    //free
    for(int iterm=0;iterm<nterms;iterm++)
      {
	nissa_free(v_o[iterm]);
	nissa_free(chi_e[iterm]);
      }
    
    nissa_free(chi_e);
    nissa_free(v_o);
    
    STOP_TIMING(quark_force_over_time);
  }
  THREADABLE_FUNCTION_END
}
