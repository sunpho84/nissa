#ifndef _SU3_OP_HPP
#define _SU3_OP_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <math.h>
#include <string.h>

#include "geometry/geometry_eo.hpp"

#include "metaprogramming/templateEnabler.hpp"

#include "routines/ios.hpp"

#include "complex.hpp"
#include "su3.hpp"

#ifdef USE_EIGEN
/// Suppress an unitialized variable warning
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
# include <Eigen/Dense>
# include <Eigen/Eigenvalues>
#pragma GCC diagnostic pop

 using namespace Eigen;

 typedef std::complex<double> scomplex_t;
 typedef scomplex_t emw_t[NCOL*NCOL];
 typedef Matrix<scomplex_t,NCOL,NCOL,RowMajor> esu3_t;
 typedef Map<esu3_t> mesu3_t;

#define SU3_ECAST(A) (mesu3_t(*(emw_t*)(A)))
#define CCAST(A) (*(scomplex_t*)(A))

#endif

namespace nissa
{
#if NCOL == 3
  //used to exponentiate for stouting
  struct hermitian_exp_ingredients
  {
    int sign;
    double c0,c1;
    double cu,su;
    double c2u,s2u;
    su3 Q,Q2;
    double u,w,theta;
    double xi0w;
    double cw;
    complex f[3];
  };
  
  CUDA_MANAGED extern int su3_sub_gr_indices[3][2];
#endif
  
  ////////////////////////////////////////////////////////////////////
  
  extern su3 gell_mann_matr[NCOL*NCOL-1];
  
  /// m=0
  template <typename U>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_put_to_zero(U&& m)
  {
    for(int ic=0;ic<NCOL;ic++)
      complex_put_to_zero(m[ic]);
  }
  
  /// m=0
  template <typename U>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_put_to_zero(U&& m)
  {
    for(int ic=0;ic<NCOL;ic++)
      color_put_to_zero(m[ic]);
  }
  
  CUDA_HOST_AND_DEVICE inline void as2t_su3_put_to_zero(as2t_su3 m) {for(size_t i=0;i<sizeof(as2t_su3)/sizeof(su3);i++) su3_put_to_zero(m[i]);}
  
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spincolor_put_to_zero(A&& a)
  {
    for(int id=0;id<NDIRAC;id++)
      color_put_to_zero(a[id]);
  }
  
  inline void colorspinspin_put_to_zero(colorspinspin m) {for(size_t ic=0;ic<NCOL;ic++) spinspin_put_to_zero(m[ic]);}
  inline void su3spinspin_put_to_zero(su3spinspin m) {for(size_t ic=0;ic<NCOL;ic++) colorspinspin_put_to_zero(m[ic]);}
  
  /// m=diag(1)
  template <typename U>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_put_to_id(U&& m)
  {
    su3_put_to_zero(m);
    
    for(int ic=0;ic<NCOL;ic++)
      m[ic][ic][RE]=1;
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_put_to_diag_color(A&& m,
			     const B& in)
  {
    su3_put_to_zero(m);
    for(int ic=0;ic<NCOL;ic++)
      complex_copy(m[ic][ic],in[ic]);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_put_to_diag_complex(A&& m,
			       const B& in)
  {
    su3_put_to_zero(m);
    for(int ic=0;ic<NCOL;ic++)
      complex_copy(m[ic][ic],in);
  }
  
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_put_to_diag_double(A&& m,
			      const double& in)
  {
    su3_put_to_zero(m);
    for(int ic=0;ic<NCOL;ic++)
      m[ic][ic][0]=in;
  }
  
  //////////////////////////////////////// Copy /////////////////////////////////////
  
  /// a=b
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_copy(A&& a,
		  const B& b)
  {
    for(int ic=0;ic<NCOL;ic++)
      complex_copy(a[ic],b[ic]);
  }
  
  /// a=b
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_copy(A&& a,
		const B& b)
  {
    for(int ic=0;ic<NCOL;ic++)
      color_copy(a[ic],b[ic]);
  }
  
  /// a=b
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void quad_su3_copy(A&& a,
		     const B& b)
  {
    for(int i=0;i<NDIM;i++)
      su3_copy(a[i],b[i]);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spincolor_copy(A&& a,
		      const B& b)
  {
    for(int i=0;i<NDIRAC;i++)
      color_copy(a[i],b[i]);
  }
  
  template <typename Out,
	    typename In>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void colorspinspin_copy(Out&& b,
			  const In& a)
  {
    for(int i=0;i<NCOL;i++)
      spinspin_copy(b[i],a[i]);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3spinspin_copy(A&& a,
			const B& b)
  {
    for(int i=0;i<NCOL;i++)
      colorspinspin_copy(a[i],b[i]);
  }
  
  //////////////////// Switch directions so to agree to ILDG ordering ////////////////////
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void quad_su3_nissa_to_ildg_reord(A&& out,
				    const B& in)
  {
    quad_su3 buff;
    quad_su3_copy(buff,in);
    
    for(int mu=0;mu<NDIM;mu++)
      su3_copy(out[(mu+NDIM-1)%NDIM],buff[mu]);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void quad_su3_ildg_to_nissa_reord(A&& out,
				    const B& in)
  {
    quad_su3 buff;
    quad_su3_copy(buff,in);
    
    for(int mu=0;mu<NDIM;mu++)
      su3_copy(out[(mu+1)%NDIM],buff[mu]);
  }
  
  ////////////////////////////////// Operations between colors //////////////////////////
  
  //just print a color
  inline void color_print(const color c)
  {
    for(size_t ic=0;ic<NCOL;ic++) printf("%+16.16lg,%+16.16lg\t",c[ic][RE],c[ic][IM]);
    master_printf("\n");
  }
  
  /// a=b+c
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_summ(A&& a,
		  const B& b,
		  const C& c)
  {
    for(int ic=0;ic<NCOL;ic++)
      complex_summ(a[ic],b[ic],c[ic]);
  }
  
  /// a=b+i*c
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_isumm(A&& a,
		  const B& b,
		  const C& c)
  {
    for(int ic=0;ic<NCOL;ic++)
      complex_isumm(a[ic],b[ic],c[ic]);
  }
  
  /// a=b-i*c
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_isubt(A&& a,
		   const B& b,
		   const C& c)
  {
    for(int ic=0;ic<NCOL;ic++)
      complex_isubt(a[ic],b[ic],c[ic]);
  }
  
  /// a=b-c
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_subt(A&& a,
		  const B& b,
		  const C& c)
  {
    for(int ic=0;ic<NCOL;ic++)
      complex_subt(a[ic],b[ic],c[ic]);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_summassign(A&& a,
			const B& b)
  {
    color_summ(a,a,b);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_subtassign(A&& a,
			const B& b)
  {
    color_subt(a,a,b);
  }
  
  CUDA_HOST_AND_DEVICE inline void color_subtassign(color a,const color b) {color_subt(a,a,b);}
  
  /// a+=i*b
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_isummassign(A&& a,
			 const B& b)
  {
    color_isumm(a,a,b);
  }
  
  /// a-=i*b
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_isubtassign(A&& a,
			 const B& b)
  {
    color_isubt(a,a,b);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_prod_double(A&& a,
			 const B& b,
			 const double& c)
  {
    for(int ic=0;ic<NCOL;ic++)
      complex_prod_double(a[ic],b[ic],c);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_prod_idouble(A&& a,
			  const B& b,
			  const C& c)
  {
    for(int ic=0;ic<NCOL;ic++)
      complex_prod_idouble(a[ic],b[ic],c);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_summ_the_prod_double(A&& a,
				  const B& b,
				  const double& c)
  {
    for(int ic=0;ic<NCOL;ic++)
      complex_summ_the_prod_double(a[ic],b[ic],c);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_summ_the_prod_idouble(A&& a,
				   const B& b,
				   const C& c)
  {
    for(int ic=0;ic<NCOL;ic++)
      complex_summ_the_prod_idouble(a[ic],b[ic],c);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_summ_the_prod_complex(A&& a,
				   const B& b,
				   const C& c)
  {
    for(int ic=0;ic<NCOL;ic++)
      complex_summ_the_prod(a[ic],b[ic],c);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_subt_the_prod_complex(A&& a,
				 const B& b,
				 const C& c)
  {
    for(size_t ic=0;ic<NCOL;ic++)
      complex_subt_the_prod(a[ic],b[ic],c);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_prod_color_delta(A&& a,
			      const B& b,
			      const size_t& c)
  {
    for(size_t ic=0;ic<NCOL;ic++)
      complex_prod_double(a[ic],b[ic],c==ic);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_linear_comb(A&& a,
			 const B& b,
			 const double& cb,
			 const C& c,
			 const double& cc)
  {
    for(size_t ic=0;ic<NCOL;ic++)
      complex_linear_comb(a[ic],b[ic],cb,c[ic],cc);
  }
  
  /// a=(b,c)
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void color_scalar_prod(A&& a,
			 const B& b,
			 const C& c)
  {
    unsafe_complex_conj1_prod(a,b[0],c[0]);
    for(int ic=1;ic<NCOL;ic++)
      complex_summ_the_conj1_prod(a,b[ic],c[ic]);
  }
  
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  double color_norm2(A&& c)
  {
    double out=complex_norm2(c[0]);
    
    for(int ic=1;ic<NCOL;ic++)
      out+=complex_norm2(c[ic]);
    
    return out;
  }
  
  //////////////////////////////////// Color and complex //////////////////////////////
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_color_prod_complex(A&& out,
			       const B& in,
			       const C& factor)
  {
    for(size_t i=0;i<NCOL;i++)
      safe_complex_prod(out[i],in[i],factor);
  }
  
  template <typename O,
	    typename I,
	    typename F>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_color_prod_complex(O&& out,
				 const I& in,
				 const F& factor)
   {
     for(int i=0;i<NCOL;i++)
      unsafe_complex_prod(out[i],in[i],factor);
  }
  
  inline void safe_color_prod_complex_conj(color out,const color in,const complex factor)
  {for(size_t i=0;i<NCOL;i++) safe_complex_conj2_prod(((complex*)out)[i],((complex*)in)[i],factor);}
  inline void unsafe_color_prod_complex_conj(color out,const color in,const complex factor)
  {for(size_t i=0;i<NCOL;i++) unsafe_complex_conj2_prod(((complex*)out)[i],((complex*)in)[i],factor);}
  
  ////////////////////////////////// Operations between su3 //////////////////////////
  
  //comput third row according to other 2
  inline void su3_build_third_row(su3 o)
  {
    CRASH_IF_NOT_3COL();
    
    unsafe_complex_conj_conj_prod(o[2][0],o[0][1],o[1][2]);
    complex_subt_the_conj_conj_prod(o[2][0],o[0][2],o[1][1]);
    unsafe_complex_conj_conj_prod(o[2][1],o[0][2],o[1][0]);
    complex_subt_the_conj_conj_prod(o[2][1],o[0][0],o[1][2]);
    unsafe_complex_conj_conj_prod(o[2][2],o[0][0],o[1][1]);
    complex_subt_the_conj_conj_prod(o[2][2],o[0][1],o[1][0]);
  }
  
  /// Just print an su3 matrix
  template <typename U>
  inline void su3_print(const U& u)
  {
    for(size_t ic1=0;ic1<NCOL;ic1++)
      {
	for(size_t ic2=0;ic2<NCOL;ic2++)
	  printf("%+16.16lg,%+16.16lg\t",u[ic1][ic2][RE],u[ic1][ic2][IM]);
	printf("\n");
      }
    printf("\n");
  }
  
  //return the trace of an su3 matrix
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_trace(A&& tr,
		 const B& m)
  {
    complex_copy(tr,m[0][0]);
    for(int ic=1;ic<NCOL;ic++)
      complex_summassign(tr,m[ic][ic]);
  }
  
  /// Returns only the real part of an su3 matrix
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  auto su3_real_trace(const A& m)
  {
    auto a=m[0][0][RE];
    
    for(int ic=1;ic<NCOL;ic++)
      a+=m[ic][ic][RE];
    
    return a;
  }
  
  //take projection of the su2 matrix over an su3 matrix
  //return the inverse modulo of the part parallel in the original matrix
  CUDA_HOST_AND_DEVICE inline double su2_part_of_su3(double &A,double &B,double &C,double &D,const su3 in,const int isub_gr)
  {
    //take indices of the subgroup
    int a=su3_sub_gr_indices[isub_gr][0];
    int b=su3_sub_gr_indices[isub_gr][1];
    
    //extract part parallel to sigmas
    A=in[a][a][RE]+in[b][b][RE];
    B=in[a][b][IM]+in[b][a][IM];
    C=in[a][b][RE]-in[b][a][RE];
    D=in[a][a][IM]-in[b][b][IM];
    
    //normalize
    double N=sqrt(A*A+B*B+C*C+D*D);
    if(fabs(N)==0) N=A=1;
    else
      {
	N=1/N;
	A*=N;
	B*=N;
	C*=N;
	D*=N;
      }
    
    return N;
  }
  
  //inverse of and su3
  CUDA_HOST_AND_DEVICE inline void su2_inv(double &a,double &b,double &c,double &d,const double a0,const double b0,const double c0,const double d0)
  {a=a0;b=-b0;c=-c0;d=-d0;}
  
  //return the overrelaxing
  CUDA_HOST_AND_DEVICE inline void su2_get_overrelaxing(double &x0,double &x1,double &x2,double &x3,const double r0,const double r1,const double r2,const double r3)
  {
    x0=2*r0*r0-1;
    x1=-2*r0*r1;
    x2=-2*r0*r2;
    x3=-2*r0*r3;
  }
  
  //multiply an su2 matrix and an su3 and assign to last
  CUDA_HOST_AND_DEVICE inline void su2_prodassign_su3(const su2 mod,const int isub_gr,su3 in)
  {
    int ic1=su3_sub_gr_indices[isub_gr][0];
    int ic2=su3_sub_gr_indices[isub_gr][1];
    
    //create the two new rows of the matrix, column by column
    for(size_t ic=0;ic<NCOL;ic++)
      {
	complex row[2];
	for(int irow=0;irow<2;irow++)
	  {
	    safe_complex_prod    (row[irow],mod[irow][0],in[ic1][ic]);
	    complex_summ_the_prod(row[irow],mod[irow][1],in[ic2][ic]);
	  }
	
	//change the two lines in the matrix
	complex_copy(in[ic1][ic],row[0]);
	complex_copy(in[ic2][ic],row[1]);
      }
  }
  
  //in the form A+B*i*sigma1+...
  CUDA_HOST_AND_DEVICE inline void su2_prodassign_su3(const double A,const double B,const double C,const double D,const int isub_gr,su3 in)
  {
    const su2 mod={{{A,D},{C,B}},{{-C,B},{A,-D}}};
    su2_prodassign_su3(mod,isub_gr,in);
  }
  
  //return the norm of U-1
  CUDA_HOST_AND_DEVICE inline double su2_nonunitarity(double A,const double B,const double C,const double D)
  {return A*A-2*A+1+B*B+C*C+D*D;}
  inline double su2_nonunitarity(su2 mod)
  {return su2_nonunitarity(mod[0][0][RE],mod[0][1][IM],mod[0][1][RE],mod[0][0][IM]);}
  
  //summ the trace to the input
  inline void su3_summ_the_trace(complex tr,const su3 m)
  {for(size_t ic=0;ic<NCOL;ic++) complex_summ(tr,tr,m[ic][ic]);}
  
  /// Return the anti-hermitian traceless part of an su3 matrix
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3_traceless_anti_hermitian_part(A&& out,
						const B& in)
  {
    double trace_im_third=0;
    for(int ic=0;ic<NCOL;ic++)
      trace_im_third+=in[ic][ic][IM];
    trace_im_third/=NCOL;
    
    for(int ic=0;ic<NCOL;ic++)
      {
	//real part of diagonal: 0
	out[ic][ic][RE]=0;
	//imag part of diagonal: subtract the trace
	out[ic][ic][IM]=in[ic][ic][IM]-trace_im_third;
	
	for(int jc=0;jc<ic;jc++)
	  {
	    //out-of-diag real part
	    out[ic][jc][RE]=-(out[jc][ic][RE]=(in[jc][ic][RE]-in[ic][jc][RE])/2);
	    //out-of-diag imag part
	    out[ic][jc][IM]=+(out[jc][ic][IM]=(in[jc][ic][IM]+in[ic][jc][IM])/2);
	  }
      }
  }
  
  //keep the antihermitian part of an su3 matrix
  CUDA_HOST_AND_DEVICE inline void su3_anti_hermitian_part(su3 out,su3 in)
  {
    for(size_t ic=0;ic<NCOL;ic++)
      {
	//real part of diagonal: 0
	out[ic][ic][RE]=0;
	//imag part of diagonal: no change
	out[ic][ic][IM]=in[ic][ic][IM];
	
	for(size_t jc=0;jc<ic;jc++)
	  {
	    //out-of-diag real part
	    out[ic][jc][RE]=-(out[jc][ic][RE]=(in[jc][ic][RE]-in[ic][jc][RE])/2);
	    //out-of-diag imag part
	    out[ic][jc][IM]=+(out[jc][ic][IM]=(in[jc][ic][IM]+in[ic][jc][IM])/2);
	  }
      }
  }
  
  //return the hermitian traceless part of an su3 matrix
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3_traceless_hermitian_part(A&& out,
					   const B& in)
  {
    double trace_re_third=0;
    for(int ic=0;ic<NCOL;ic++) trace_re_third+=in[ic][ic][RE];
    trace_re_third/=NCOL;
    
    for(int ic=0;ic<NCOL;ic++)
      {
	//imag part of diagonal: 0
	out[ic][ic][IM]=0;
	//real part of diagonal: subtract the trace
	out[ic][ic][RE]=in[ic][ic][RE]-trace_re_third;
	
	for(int jc=0;jc<ic;jc++)
	  {
	    //out-of-diag real part
	    out[ic][jc][RE]=+(out[jc][ic][RE]=(in[jc][ic][RE]+in[ic][jc][RE])/2);
	    //out-of-diag imag part
	    out[ic][jc][IM]=-(out[jc][ic][IM]=(in[jc][ic][IM]-in[ic][jc][IM])/2);
	  }
      }
  }
  
  //calculate the determinant of an su3 matrix
  CUDA_HOST_AND_DEVICE inline void su3_det(complex d,const su3 U)
  {
#if NCOL == 3
    complex a;
    
    unsafe_complex_prod(a,U[1][1],U[2][2]);
    complex_subt_the_prod(a,U[1][2],U[2][1]);
    unsafe_complex_prod(d,U[0][0],a);
    
    unsafe_complex_prod(a,U[1][2],U[2][0]);
    complex_subt_the_prod(a,U[1][0],U[2][2]);
    complex_summ_the_prod(d,U[0][1],a);
    
    unsafe_complex_prod(a,U[1][0],U[2][1]);
    complex_subt_the_prod(a,U[1][1],U[2][0]);
    complex_summ_the_prod(d,U[0][2],a);
#else
    matrix_determinant(d,(complex*)U,NCOL);
#endif
  }
  
  //calculate the real part of the determinant of an su3 matrix
  CUDA_HOST_AND_DEVICE inline double su3_real_det(const su3 u)
  {
#if NCOL == 3
    return
      u[0][2][IM]*(u[1][1][RE]*u[2][0][IM]+u[1][1][IM]*u[2][0][RE]-u[1][0][RE]*u[2][1][IM]-u[1][0][IM]*u[2][1][RE])+
      u[0][2][RE]*(u[1][1][IM]*u[2][0][IM]-u[1][1][RE]*u[2][0][RE]-u[1][0][IM]*u[2][1][IM]+u[1][0][RE]*u[2][1][RE])+
      u[0][1][IM]*(u[1][0][RE]*u[2][2][IM]-u[1][2][RE]*u[2][0][IM]-u[1][2][IM]*u[2][0][RE]+u[1][0][IM]*u[2][2][RE])+
      u[0][1][RE]*(u[1][0][IM]*u[2][2][IM]-u[1][2][IM]*u[2][0][IM]+u[1][2][RE]*u[2][0][RE]-u[1][0][RE]*u[2][2][RE])+
      u[0][0][IM]*(u[1][2][RE]*u[2][1][IM]+u[1][2][IM]*u[2][1][RE]-u[1][1][RE]*u[2][2][IM]-u[1][1][IM]*u[2][2][RE])+
      u[0][0][RE]*(u[1][2][IM]*u[2][1][IM]-u[1][2][RE]*u[2][1][RE]-u[1][1][IM]*u[2][2][IM]+u[1][1][RE]*u[2][2][RE]);
#else
    complex d;
    su3_det(d,u);
    return d[RE];
#endif
  }
  
  /// Return the hemitian su3 matrix
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3_hermitian(A&& out,const B& in)
  {
    for(int ic_in=0;ic_in<NCOL;ic_in++)
      for(int ic_out=0;ic_out<NCOL;ic_out++)
	complex_conj(out[ic_in][ic_out],in[ic_out][ic_in]);
  }
  
  inline void safe_su3_hermitian(su3 out,const su3 in)
  {
    su3 tmp;
    unsafe_su3_hermitian(tmp,in);
    su3_copy(out,tmp);
  }
  
  //return the transposed su3 matrix
  inline void unsafe_su3_transpose(su3 out,const su3 in)
  {
    for(size_t ic_in=0;ic_in<NCOL;ic_in++)
      for(size_t ic_out=0;ic_out<NCOL;ic_out++)
	complex_copy(out[ic_in][ic_out],in[ic_out][ic_in]);
  }
  inline void safe_su3_transpose(su3 out,const su3 in)
  {
    su3 tmp;
    unsafe_su3_transpose(tmp,in);
    su3_copy(out,tmp);
  }
  
  /// a=b+c
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_summ(A&& a,
		const B& b,
		const C& c)
  {
    for(int ic=0;ic<NCOL;ic++)
      color_summ(a[ic],b[ic],c[ic]);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3_summ_su3_dag(A&& a,
			       const B& b,
			       const C& c)
  {
    for(int i=0;i<NCOL;i++)
      for(int j=0;j<NCOL;j++)
	complex_summ_conj2(a[i][j],b[i][j],c[j][i]);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_summassign_su3_dag(A&& a,
			      const B& b)
  {
    unsafe_su3_summ_su3_dag(a,a,b);
  }
  
  /// a+=b
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_summassign(A&& a,
		      const B& b)
  {
    su3_summ(a,a,b);
  }
  
  inline void su3_summ_real(su3 a,const su3 b,const double c)
  {su3_copy(a,b);for(size_t i=0;i<NCOL;i++) a[i][i][0]=b[i][i][0]+c;}
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_subt(A&& a,
		const B& b,
		const C& c)
  {
    for(int ic=0;ic<NCOL;ic++)
      color_subt(a[ic],b[ic],c[ic]);
  }
  
  /// a+=b
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_subtassign(A&& a,
		      const B& b)
  {
    su3_subt(a,a,b);
  }
  
  inline void su3_subt_complex(su3 a,const su3 b,const complex c) {su3_copy(a,b);for(size_t i=0;i<NCOL;i++) complex_subt(a[i][i],b[i][i],c);}
  CUDA_HOST_AND_DEVICE inline void unsafe_su3_subt_su3_dag(su3 a,const su3 b,const su3 c)
  {
    for(size_t i=0;i<NCOL;i++)
      for(size_t j=0;j<NCOL;j++)
	complex_subt_conj2(a[i][j],b[i][j],c[j][i]);
  }
  inline void su3_subtassign_su3_dag(su3 a,const su3 b) {unsafe_su3_subt_su3_dag(a,a,b);}
  
  /// Product of two su3 matrixes
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3_prod_su3(A&& a,
			   const B& b,
			   const C& c,
			   const size_t nr_max=NCOL)
  {
    for(size_t ir_out=0;ir_out<nr_max;ir_out++)
      for(size_t ic_out=0;ic_out<NCOL;ic_out++)
	{
	  unsafe_complex_prod(a[ir_out][ic_out],b[ir_out][0],c[0][ic_out]);
	  for(int itemp=1;itemp<NCOL;itemp++)
	    complex_summ_the_prod(a[ir_out][ic_out],b[ir_out][itemp],c[itemp][ic_out]);
	}
  }

  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_su3_prod_su3(A&& a,
			 const B& b,
			 const C& c)
  {
    su3 d;
    unsafe_su3_prod_su3(d,b,c);
    su3_copy(a,d);
  }
  
  inline void su3_prodassign_su3(su3 a,const su3 b) {safe_su3_prod_su3(a,a,b);}
  
  /// a+=b*c
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_summ_the_prod_su3(A&& a,
			     const B& b,
			     const C& c)
  {
    for(int ir_out=0;ir_out<NCOL;ir_out++)
      for(int ic_out=0;ic_out<NCOL;ic_out++)
	for(int itemp=0;itemp<NCOL;itemp++)
	  complex_summ_the_prod(a[ir_out][ic_out],b[ir_out][itemp],c[itemp][ic_out]);
  }
  
  /// a+=b*c^+
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_summ_the_prod_su3_dag(A&& a,
				 const B& b,
				 const C& c)
  {
    for(int ir_out=0;ir_out<NCOL;ir_out++)
      for(int ic_out=0;ic_out<NCOL;ic_out++)
	for(int itemp=0;itemp<NCOL;itemp++)
	  complex_summ_the_conj2_prod(a[ir_out][ic_out],b[ir_out][itemp],c[ic_out][itemp]);
  }
  
  /// a+=b^*c
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_summ_the_dag_prod_su3(A&& a,
				 const B& b,
				 const C& c)
  {
    for(int ir_out=0;ir_out<NCOL;ir_out++)
      for(int ic_out=0;ic_out<NCOL;ic_out++)
	for(int itemp=0;itemp<NCOL;itemp++)
	  complex_summ_the_conj1_prod(a[ir_out][ic_out],b[itemp][ir_out],c[itemp][ic_out]);
  }
  
  /// Product of two su3 matrixes
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3_dag_prod_su3(A&& a,
			       const B& b,
			       const C& c,
			       const int nr_max=NCOL)
  {
    for(int ir_out=0;ir_out<nr_max;ir_out++)
      for(int ic_out=0;ic_out<NCOL;ic_out++)
	{
	  unsafe_complex_conj1_prod(a[ir_out][ic_out],b[0][ir_out],c[0][ic_out]);
	  for(int itemp=1;itemp<NCOL;itemp++)
	    complex_summ_the_conj1_prod(a[ir_out][ic_out],b[itemp][ir_out],c[itemp][ic_out]);
	}
  }
  
  inline void safe_su3_dag_prod_su3(su3 a,const su3 b,const su3 c) {su3 d;unsafe_su3_dag_prod_su3(d,b,c);su3_copy(a,d);}
  inline void su3_dag_summ_the_prod_su3(su3 a,const su3 b,const su3 c)
  {
    for(size_t ir_out=0;ir_out<NCOL;ir_out++)
      for(size_t ic_out=0;ic_out<NCOL;ic_out++)
	for(size_t itemp=0;itemp<NCOL;itemp++)
	  complex_summ_the_conj1_prod(a[ir_out][ic_out],b[itemp][ir_out],c[itemp][ic_out]);
  }
  
  /// Product of two su3 matrixes
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3_prod_su3_dag(A&& a,
			       const B& b,
			       const C& c,
			       const int nr_max=NCOL)
  {
    for(int ir_out=0;ir_out<nr_max;ir_out++)
      for(int ic_out=0;ic_out<NCOL;ic_out++)
	{
	  unsafe_complex_conj2_prod(a[ir_out][ic_out],b[ir_out][0],c[ic_out][0]);
	  for(int jc=1;jc<NCOL;jc++)
	    complex_summ_the_conj2_prod(a[ir_out][ic_out],b[ir_out][jc],c[ic_out][jc]);
	}
  }
  
  CUDA_HOST_AND_DEVICE inline void safe_su3_prod_su3_dag(su3 a,const su3 b,const su3 c) {su3 d;unsafe_su3_prod_su3_dag(d,b,c);su3_copy(a,d);}
  
  /// subtract the product
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_subt_the_prod_su3_dag(A&& a,const B& b,const C& c)
  {
    for(size_t ir_out=0;ir_out<NCOL;ir_out++)
      for(size_t ic_out=0;ic_out<NCOL;ic_out++)
	for(size_t jc=0;jc<NCOL;jc++)
	  complex_subt_the_conj2_prod(a[ir_out][ic_out],b[ir_out][jc],c[ic_out][jc]);
  }
  
  /// Trace of the product of two su3 matrices
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  double real_part_of_trace_su3_prod_su3_dag(A&& a,
					     const B& b)
  {
    double t=0.0;
    
    for(int ic1=0;ic1<NCOL;ic1++)
      for(int ic2=0;ic2<NCOL;ic2++)
	t+=real_part_of_complex_scalar_prod(a[ic1][ic2],b[ic1][ic2]);
    
    return t;
  }
  
  //Trace of the product of two su3 matrices
  CUDA_HOST_AND_DEVICE inline void trace_su3_prod_su3(complex t,const su3 a,const su3 b)
  {
    complex_put_to_zero(t);
    for(size_t ic1=0;ic1<NCOL;ic1++)
      for(size_t ic2=0;ic2<NCOL;ic2++)
	complex_summ_the_prod(t,a[ic1][ic2],b[ic2][ic1]);
  }
  
  /// Product of two su3 matrixes
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3_dag_prod_su3_dag(A&& a,
				   const B& b,
				   const C& c,
				   const int nr_max=NCOL)
  {
    for(int ir_out=0;ir_out<nr_max;ir_out++)
      for(int ic_out=0;ic_out<NCOL;ic_out++)
	{
	  unsafe_complex_conj_conj_prod(a[ir_out][ic_out],b[0][ir_out],c[ic_out][0]);
	  for(int itemp=1;itemp<NCOL;itemp++)
	    complex_summ_the_conj_conj_prod(a[ir_out][ic_out],b[itemp][ir_out],c[ic_out][itemp]);
	}
  }
  
  inline void safe_su3_dag_prod_su3_dag(su3 a,const su3 b,const su3 c)
  {su3 d;unsafe_su3_dag_prod_su3_dag(d,b,c);su3_copy(a,d);}
  
  /// Product of an su3 matrix by a complex
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3_prod_complex(A&& a,
			       const B& b,
			       const C& c)
  {
    for(int ic=0;ic<NCOL;ic++)
      unsafe_color_prod_complex(a[ic],b[ic],c);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_su3_prod_complex(A&& a,
			       const B& b,
			       const C& c)
  {
    su3 d;
    unsafe_su3_prod_complex(d,b,c);
    su3_copy(a,d);
  }
  
  //product of an su3 matrix by a complex
  inline void unsafe_su3_prod_complex_conj(su3 a,const su3 b,const complex c)
  {for(size_t ic=0;ic<NCOL;ic++) unsafe_color_prod_complex_conj(a[ic],b[ic],c);}
  inline void safe_su3_prod_complex_conj(su3 a,const su3 b,const complex c)
  {su3 d;unsafe_su3_prod_complex_conj(d,b,c);su3_copy(a,d);}
  
  /// Product of an su3 matrix by a real
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_prod_double(A&& a,
		       const B& b,
		       const double& r)
  {
    for(int ic=0;ic<NCOL;ic++)
      color_prod_double(a[ic],b[ic],r);
  }
  
  template <typename A>
  CUDA_HOST_AND_DEVICE inline void su3_prodassign_double(A&& a,
							 const double& r)
  {
    su3_prod_double(a,a,r);
  }
  
  /// Hermitian of su3 matrix times a real
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3_hermitian_prod_double(A&& a,
					const B& b,
					const double& r)
  {
    for(int i=0;i<NCOL;i++)
      for(int j=0;j<NCOL;j++)
	{
	  a[i][j][RE]=+r*b[j][i][RE];
	  a[i][j][IM]=-r*b[j][i][IM];
	}
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_su3_hermitian_prod_double(A&& a,
				      const B& b,
				      const double& r)
  {
    for(int i=0;i<NCOL;i++)
      {
	a[i][i][RE]=+r*b[i][i][RE];
	a[i][i][IM]=-r*b[i][i][IM];
	
	for(int j=i+1;j<NCOL;j++)
	  {
	    const double a_i_j_RE=+r*b[j][i][RE];
	    const double a_i_j_IM=-r*b[j][i][IM];
	    a[j][i][RE]=+r*b[i][j][RE];
	    a[j][i][IM]=-r*b[i][j][IM];
	    a[i][j][RE]=a_i_j_RE;
	    a[i][j][IM]=a_i_j_IM;
	  }
      }
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_summ_the_hermitian_prod_double(A&& a,
					  const B& b,
					  const double& r)
  {
    for(int i=0;i<NCOL;i++)
      for(int j=0;j<NCOL;j++)
	{
	  a[i][j][RE]+= r*b[j][i][RE];
	  a[i][j][IM]+=-r*b[j][i][IM];
	}
  }
  
  /// Summ the prod of su3 with imag
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_prod_idouble(A&& a,
			const B& b,
			const double& r)
  {
    for(int i=0;i<NCOL;i++)
      color_prod_idouble(a[i],b[i],r);
  }
  
  inline void su3_prodassign_idouble(su3 a,const double r) {su3_prod_idouble(a,a,r);}
  
  /// a+=b*i*r
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_summ_the_prod_idouble(A&& a,
				 const B& b,
				 const double& r)
  {
    for(int i=0;i<NCOL;i++)
      color_summ_the_prod_idouble(a[i],b[i],r);
  }
  
  //summ the prod of su3 with real
  CUDA_HOST_AND_DEVICE inline void su3_summ_the_prod_double(su3 a,const su3 b,const double r) {for(size_t i=0;i<NCOL;i++) color_summ_the_prod_double(a[i],b[i],r);}
  
  //summ the prod of the dag su3 with real
  inline void su3_dag_summ_the_prod_double(su3 a,const su3 b,const double r)
  {
    for(size_t i=0;i<NCOL;i++)
      for(size_t j=0;j<NCOL;j++)
	{
	  a[i][j][RE]+=b[j][i][RE]*r;
	  a[i][j][IM]-=b[j][i][IM]*r;
	}
  }
  
  //combine linearly two su3 elements
  inline void su3_linear_comb(su3 a,const su3 b,const double cb,su3 c,const double cc) {for(size_t ic=0;ic<NCOL;ic++) color_linear_comb(a[ic],b[ic],cb,c[ic],cc);}
  
  /// Summ the prod of su3 with complex
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_summ_the_prod_complex(A&& a,
				 const B& b,
				 const C& c)
  {
    for(int ic=0;ic<NCOL;ic++)
      color_summ_the_prod_complex(a[ic],b[ic],c);
  }
  
  //calculate explicitely the inverse
  CUDA_HOST_AND_DEVICE inline void unsafe_su3_explicit_inverse(su3 invU,const su3 U)
  {
    CRASH_IF_NOT_3COL();
    
    complex det,rec_det;
    su3_det(det,U);
    complex_reciprocal(rec_det,det);
    
    unsafe_complex_prod(invU[0][0],U[1][1],U[2][2]);
    unsafe_complex_prod(invU[1][0],U[1][2],U[2][0]);
    unsafe_complex_prod(invU[2][0],U[1][0],U[2][1]);
    
    unsafe_complex_prod(invU[0][1],U[0][2],U[2][1]);
    unsafe_complex_prod(invU[1][1],U[0][0],U[2][2]);
    unsafe_complex_prod(invU[2][1],U[0][1],U[2][0]);
    
    unsafe_complex_prod(invU[0][2],U[0][1],U[1][2]);
    unsafe_complex_prod(invU[1][2],U[0][2],U[1][0]);
    unsafe_complex_prod(invU[2][2],U[0][0],U[1][1]);
    
    
    complex_subt_the_prod(invU[0][0],U[1][2],U[2][1]);
    complex_subt_the_prod(invU[1][0],U[1][0],U[2][2]);
    complex_subt_the_prod(invU[2][0],U[1][1],U[2][0]);
    
    complex_subt_the_prod(invU[0][1],U[0][1],U[2][2]);
    complex_subt_the_prod(invU[1][1],U[0][2],U[2][0]);
    complex_subt_the_prod(invU[2][1],U[0][0],U[2][1]);
    
    complex_subt_the_prod(invU[0][2],U[0][2],U[1][1]);
    complex_subt_the_prod(invU[1][2],U[0][0],U[1][2]);
    complex_subt_the_prod(invU[2][2],U[0][1],U[1][0]);
    
    for(size_t icol1=0;icol1<NCOL;icol1++)
      for(size_t icol2=0;icol2<NCOL;icol2++)
	safe_complex_prod(invU[icol1][icol2],invU[icol1][icol2],rec_det);
  }
  inline void safe_su3_explicit_inverse(su3 invU,const su3 U)
  {su3 tempU;unsafe_su3_explicit_inverse(tempU,U);su3_copy(invU,tempU);}
  
  //summ of the squared norm of the entries
  CUDA_HOST_AND_DEVICE inline double su3_norm2(const su3 U)
  {
    double norm2=0;
    for(size_t ic=0;ic<NCOL;ic++) norm2+=color_norm2(U[ic]);
    
    return norm2;
  }
  
  //compute the square root of y numerically
  inline void unsafe_su3_sqrt(su3 x,const su3 y)
  {
    su3_put_to_id(x);
    
    double err;
    do
      {
	//update
	su3 t;
	unsafe_su3_explicit_inverse(t,x);
	safe_su3_prod_su3(t,y,t);
	su3_subt(t,t,x);
	su3_summ_the_prod_double(x,t,0.5);
	
	//compute the error
	unsafe_su3_prod_su3(t,x,x);
	su3_subt(t,y,t);
	err=sqrt(su3_norm2(t));
      }
    while(err>1.e-15);
  }
  inline void safe_su3_sqrt(su3 x,const su3 y) {su3 t;unsafe_su3_sqrt(t,y);su3_copy(x,t);}
  
  //exponentiate an su3 matrix through taylor expansion
  inline void unsafe_su3_taylor_exponentiate(su3 out,const su3 in,size_t order)
  {
    //1st terms
    double coef=1;
    su3 temp;
    su3_copy(temp,in);
    
    //order 0+1
    su3_put_to_id(out);
    su3_summassign(out,in);
    
    //higher orders
    for(size_t iorder=2;iorder<=order;iorder++)
      {
	safe_su3_prod_su3(temp,temp,in);
	coef/=iorder;
	su3_summ_the_prod_double(out,temp,coef);
      }
  }
  
  //unitarize by orthonormalizing the rows
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_unitarize_orthonormalizing(A&& o,
				      const B& i)
  {
    CRASH_IF_NOT_3COL();
    
    //compute the squared norm of row 0
    const double row0_norm2=color_norm2(i[0]);
    
    //compute the scalar product between row 1 and 0
    complex row10_sc_prod;
    unsafe_complex_conj2_prod(  row10_sc_prod,i[1][0],i[0][0]);
    complex_summ_the_conj2_prod(row10_sc_prod,i[1][1],i[0][1]);
    complex_summ_the_conj2_prod(row10_sc_prod,i[1][2],i[0][2]);
    
    //orthonormalize row 1
    complex f;
    complex_prod_double(f,row10_sc_prod,1/row0_norm2);
    
    for(size_t c=0;c<NCOL;c++)
      {
	for(size_t ri=0;ri<2;ri++) o[1][c][ri]=i[1][c][ri];
	complex_subt_the_prod(o[1][c],f,i[0][c]);
      }
    
    const double row0_norm=1/sqrt(row0_norm2);
    const double row1_norm=1/sqrt(complex_norm2(o[1][0])+complex_norm2(o[1][1])+complex_norm2(o[1][2]));
    
    //normalize the rows
    for(size_t c=0;c<NCOL;c++)
      for(size_t ri=0;ri<2;ri++)
	{
	  o[0][c][ri]=row0_norm*i[0][c][ri];
	  o[1][c][ri]=row1_norm*o[1][c][ri];
	}
    
    //orthonormalize the third row
    unsafe_complex_conj_conj_prod(o[2][0],o[0][1],o[1][2]);
    complex_subt_the_conj_conj_prod(o[2][0],o[0][2],o[1][1]);
    unsafe_complex_conj_conj_prod(o[2][1],o[0][2],o[1][0]);
    complex_subt_the_conj_conj_prod(o[2][1],o[0][0],o[1][2]);
    unsafe_complex_conj_conj_prod(o[2][2],o[0][0],o[1][1]);
    complex_subt_the_conj_conj_prod(o[2][2],o[0][1],o[1][0]);
  }
  
  void su3_find_overrelaxed(su3 out,const su3 in,const su3 staple,int nov_hits);
  void su3_overrelax(su3 out,const su3 in,const double w,const double *coeff,const int ord);
  
  //overrelax the link using approximated exponentiation
  inline void su3_overrelax(su3 out,const su3 in,const double w)
  {
    const double coeff[5]={1,w,w*(w-1)/2,w*(w-1)*(w-2)/6,w*(w-1)*(w-2)*(w-3)/24};
    su3_overrelax(out,in,w,coeff,5);
  }
  
  /// Unitarize an su3 matrix by taking explicitely the inverse and averaging with it
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_unitarize_explicitly_inverting(su3& new_link,
					  const su3& prop_link,
					  const double tol=1e-15)
  {
    su3 inv,temp_link;
    double gamma,residue;
    
    su3_copy(temp_link,prop_link);
    
    do
      {
	unsafe_su3_explicit_inverse(inv,temp_link);
	gamma=sqrt(su3_norm2(inv)/su3_norm2(temp_link));
	
	//average U and U^-1^+
	residue=0;
	for(size_t icol1=0;icol1<NCOL;icol1++)
	  for(size_t icol2=0;icol2<NCOL;icol2++)
	    {
	      new_link[icol1][icol2][RE]=0.5*(temp_link[icol1][icol2][RE]*gamma+inv[icol2][icol1][RE]/gamma);
	      new_link[icol1][icol2][IM]=0.5*(temp_link[icol1][icol2][IM]*gamma-inv[icol2][icol1][IM]/gamma);
	      for(size_t ri=0;ri<2;ri++)
		{
		  double diff=new_link[icol1][icol2][ri]-temp_link[icol1][icol2][ri];
		  residue+=diff*diff;
		}
	    }
	
	su3_copy(temp_link,new_link);
	
	residue=sqrt(residue);
      }
    while(residue>tol);
    
    //divide by third root of det
    complex det,fact;
    su3_det(det,new_link);
    complex_pow(fact,det,-1.0/NCOL);
    safe_su3_prod_complex(new_link,new_link,fact);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_unitarize_explicitly_inverting(A&& _new_link,
					  const B& _prop_link,
					  const double tol=1e-15)
  {
    su3 prop_link;
    su3_copy(prop_link,_prop_link);
    
    su3 new_link;
    su3_unitarize_explicitly_inverting(new_link,prop_link,tol);
    su3_copy(_new_link,new_link);
  }
  
  //perform an iteration of maximal projection trace
  //by maximing the trace over three different subgroups of su3
  inline void su3_unitarize_maximal_trace_projecting_iteration_slow(su3 U,const su3 M)
  {
    //loop over the three subgroups
    for(size_t isub_gr=0;isub_gr<NCOL;isub_gr++)
      {
	//compute the product
	su3 prod;
	unsafe_su3_prod_su3_dag(prod,U,M);
	
	//take the subgroup isub_gr
	double a,b,c,d;
	su2_part_of_su3(a,b,c,d,prod,isub_gr);
	
	//invert
	su2_inv(a,b,c,d, a,b,c,d);
	
	//modify the subgroup
	su2_prodassign_su3(a,b,c,d,isub_gr,U);
      }
  }
  
  //perform an iteration of maximal projection trace
  //by maximing the trace over three different subgroups of su3
  CUDA_HOST_AND_DEVICE inline double su3_unitarize_maximal_trace_projecting_iteration(su3 U,const su3 M)
  {
    //compute the product
    su3 prod;
    unsafe_su3_prod_su3_dag(prod,U,M);
    
    //master_printf("isub -1, %16.16lg\n",su3_real_trace(prod));
    //loop over the three subgroups
    for(size_t isub_gr=0;isub_gr<NCOL;isub_gr++)
      {
	//take the subgroup isub_gr
	double a,b,c,d;
	su2_part_of_su3(a,b,c,d,prod,isub_gr);
	
	//take the inverse
	su2_inv(a,b,c,d,a,b,c,d);
	
	//modify the subgroup
	su2_prodassign_su3(a,b,c,d,isub_gr,U);
	su2_prodassign_su3(a,b,c,d,isub_gr,prod);
	
	//master_printf("isub %d, %16.16lg\n",isub_gr,su3_real_trace(prod));
      }
    
    return su3_real_trace(prod);
  }
  
  //perform maximal projection trace up to reaching the machine precision
  CUDA_HOST_AND_DEVICE void su3_unitarize_maximal_trace_projecting(su3 out,const su3 M,const double precision=5e-15,int niter_max=20000);
  CUDA_HOST_AND_DEVICE inline void su3_unitarize_maximal_trace_projecting(su3 out,const double precision=5e-15,int niter_max=20000)
  {su3_unitarize_maximal_trace_projecting(out,out,precision,niter_max);}
  
  void su3_find_cooled_eo_conf(su3 u,eo_ptr<quad_su3> eo_conf,int par,int ieo,int mu);
  void su3_find_cooled_lx_conf(su3 u,quad_su3 *lx_conf,int ivol,int mu);
  
  ////////////////////// products between su3 and color //////////////////
  
  //product of an su3 matrix by a color vector
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3_prod_color(A&& a,
			     const B& b,
			     const C& c)
  {
    for(int c1=0;c1<NCOL;c1++)
      {
	unsafe_complex_prod(a[c1],b[c1][0],c[0]);
	for(int c2=1;c2<NCOL;c2++)
	  complex_summ_the_prod(a[c1],b[c1][c2],c[c2]);
      }
  }
  
  void su3_unitarize_with_sqrt(su3 out,const su3 in);
  
  //product of an su3 matrix by a color vector
  CUDA_HOST_AND_DEVICE inline void unsafe_single_su3_prod_single_color(single_color a,const single_su3 b,const single_color c)
  {
    for(size_t c1=0;c1<NCOL;c1++)
      {
	unsafe_single_complex_prod(a[c1],b[c1][0],c[0]);
	for(size_t c2=1;c2<NCOL;c2++) single_complex_summ_the_prod(a[c1],b[c1][c2],c[c2]);
      }
  }
  
  //safe prod
  CUDA_HOST_AND_DEVICE inline void safe_su3_prod_color(color a,const su3 b,const color c) {color t;unsafe_su3_prod_color(t,b,c);color_copy(a,t);}
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_summ_the_prod_color(A&& a,
			       const B& b,
			       const C& c)
  {
    for(int c1=0;c1<NCOL;c1++)
      for(int c2=0;c2<NCOL;c2++)
	complex_summ_the_prod(a[c1],b[c1][c2],c[c2]);
  }
  
  CUDA_HOST_AND_DEVICE inline void single_su3_summ_the_prod_single_color(single_color a,const single_su3 b,const single_color c)
  {for(size_t c1=0;c1<NCOL;c1++) for(size_t c2=0;c2<NCOL;c2++) single_complex_summ_the_prod(a[c1],b[c1][c2],c[c2]);}
  
  //subt
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_subt_the_prod_color(A&& a,
			       const B& b,
			       const C& c)
  {
    for(int c1=0;c1<NCOL;c1++)
      for(int c2=0;c2<NCOL;c2++)
	complex_subt_the_prod(a[c1],b[c1][c2],c[c2]);
  }
  
  // a=b^*c
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3_dag_prod_color(A&& a,
				 const B& b,
				 const C& c)
  {
    for(int c1=0;c1<NCOL;c1++)
      {
	unsafe_complex_conj1_prod(a[c1],b[0][c1],c[0]);
	for(int c2=1;c2<NCOL;c2++)
	  complex_summ_the_conj1_prod(a[c1],b[c2][c1],c[c2]);
      }
  }
  
  //safe dag prod
  inline void safe_su3_dag_prod_color(color a,const su3 b,const color c)
  {color t;unsafe_su3_dag_prod_color(t,b,c);color_copy(a,t);}
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_dag_summ_the_prod_color(A&& a,
				   const B& b,
				   const C& c)
  {
    for(int c1=0;c1<NCOL;c1++)
      for(int c2=0;c2<NCOL;c2++)
	complex_summ_the_conj1_prod(a[c1],b[c2][c1],c[c2]);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void single_su3_dag_summ_the_prod_single_color(A&& a,
						 const B& b,
						 const C& c)
  {
    for(size_t c1=0;c1<NCOL;c1++)
      for(size_t c2=0;c2<NCOL;c2++)
	single_complex_summ_the_conj1_prod(a[c1],b[c2][c1],c[c2]);
  }
  
  //subt dag
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_dag_subt_the_prod_color(A&& a,
				   const B& b,
				   const C& c)
  {
    for(int c1=0;c1<NCOL;c1++)
      for(int c2=0;c2<NCOL;c2++)
	complex_subt_the_conj1_prod(a[c1],b[c2][c1],c[c2]);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void single_su3_dag_subt_the_prod_single_color(A&& a,
						 const B& b,
						 const C& c)
  {
    for(size_t c1=0;c1<NCOL;c1++)
      for(size_t c2=0;c2<NCOL;c2++)
	single_complex_subt_the_conj1_prod(a[c1],b[c2][c1],c[c2]);
  }
  
  //////////////////////////////////////////// color prod su3 ///////////////////////////////////////////
  
  //prod 
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_color_prod_su3(A&& a,
			     const B& b,
			     const C& c)
  {
    for(size_t c1=0;c1<NCOL;c1++)
      {
	unsafe_complex_prod(a[c1],b[0],c[0][c1]);
	for(size_t c2=1;c2<NCOL;c2++) complex_summ_the_prod(a[c1],b[c2],c[c2][c1]);
      }
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_color_prod_su3(A&& a,
			   const B& b,
			   const C& c)
  {
    color t;
    unsafe_color_prod_su3(t,b,c);
    color_copy(a,t);
  }
  
  //dag
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_color_prod_su3_dag(A&& a,
				 const B& b,
				 const C& c)
  {
    for(size_t c1=0;c1<NCOL;c1++)
      {
	unsafe_complex_conj2_prod(a[c1],b[0],c[c1][0]);
	for(size_t c2=1;c2<NCOL;c2++)
	  complex_summ_the_conj2_prod(a[c1],b[c2],c[c1][c2]);
      }
  }
  
  //////////////////////////////// Operations between halfspincolors /////////////////////
  
  CUDA_HOST_AND_DEVICE inline double halfspincolor_scal_prod(halfspincolor a,const halfspincolor b)
  {
    double out=0;
    for(int id=0;id<NDIRAC/2;id++)
      for(int icol=0;icol<NCOL;icol++)
	for(int ri=0;ri<2;ri++)
	  out+=a[id][icol][ri]*b[id][icol][ri];
    
    return out;
  }
  CUDA_HOST_AND_DEVICE inline double halfspincolor_norm2(halfspincolor a)
  {
    return halfspincolor_scal_prod(a,a);
  }
  
  CUDA_HOST_AND_DEVICE inline void halfspincolor_summ_the_prod_double(halfspincolor a,const halfspincolor b,const halfspincolor c,const double d)
  {
    for(int id=0;id<NDIRAC/2;id++)
      for(int icol=0;icol<NCOL;icol++)
	for(int ri=0;ri<2;ri++)
	  a[id][icol][ri]=b[id][icol][ri]+c[id][icol][ri]*d;
  }
  CUDA_HOST_AND_DEVICE inline void halfspincolor_summ_the_prod_double(halfspincolor a,const halfspincolor b,const double c)
  {halfspincolor_summ_the_prod_double(a,a,b,c);}
  
  CUDA_HOST_AND_DEVICE inline void halfspincolor_copy(halfspincolor a,const halfspincolor b)
  {memcpy(a,b,sizeof(halfspincolor));}
  
  CUDA_HOST_AND_DEVICE inline void halfspincolor_put_to_zero(halfspincolor a)
  {memset(a,0,sizeof(halfspincolor));}
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_halfspincolor_halfspincolor_times_halfspincolor(A&& a,
							      const B& b,
							      const C& c,
							      const int& offset)
  {
    for(int id_out=0;id_out<NDIRAC/2;id_out++)
      for(int ic_out=0;ic_out<NCOL;ic_out++)
	{
	  complex_put_to_zero(a[id_out+offset][ic_out]);
	  for(int id_in=0;id_in<NDIRAC/2;id_in++)
	    for(int ic_in=0;ic_in<NCOL;ic_in++)
	      complex_summ_the_prod(a[id_out+offset][ic_out],b[id_out][ic_out][id_in][ic_in],c[id_in+offset][ic_in]);
	}
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_halfspincolor_halfspincolor_dag_times_halfspincolor(A&& a,
								  const B& b,
								  const C& c,
								  const int& offset)
  {
    for(int id_out=0;id_out<NDIRAC/2;id_out++)
      for(int ic_out=0;ic_out<NCOL;ic_out++)
	{
	  complex_put_to_zero(a[id_out+offset][ic_out]);
	  for(int id_in=0;id_in<NDIRAC/2;id_in++)
	    for(int ic_in=0;ic_in<NCOL;ic_in++)
	      complex_summ_the_conj1_prod(a[id_out+offset][ic_out],b[id_in][ic_in][id_out][ic_out],c[id_in+offset][ic_in]);
	}
  }
  
  ////////////////////////////////// Operations between spincolor ////////////////////////
  
  //just print a spincolor
  inline void spincolor_print(spincolor c)
  {
    for(size_t id=0;id<NDIRAC;id++)
      {
	for(size_t ic=0;ic<NCOL;ic++) printf("%+16.16lg,%+16.16lg\t",c[id][ic][0],c[id][ic][1]);
	master_printf("\n");
      }
    master_printf("\n");
  }
  
  //summ two spincolors
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spincolor_summ(A&& a,
		      const B& b,
		      const C& c)
  {
    for(int i=0;i<NDIRAC;i++)
      color_summ(a[i],b[i],c[i]);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spincolor_summassign(A&& a,
			    const B& b)
  {
    spincolor_summ(a,a,b);
  }
  
  //subtract two spincolors
  CUDA_HOST_AND_DEVICE inline void spincolor_subt(spincolor a,const spincolor b,const spincolor c) {for(size_t i=0;i<NDIRAC;i++) color_subt(a[i],b[i],c[i]);}
  inline void spincolor_subtassign(spincolor a,const spincolor b) {spincolor_subt(a,a,b);}
  
  //spincolor*real
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spincolor_prod_double(A&& a,
			     const B& b,
			     const double& factor)
  {
    for(int i=0;i<NDIRAC;i++)
      color_prod_double(a[i],b[i],factor);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spincolor_summ_the_prod_double(A&& a,
				      const B& b,
				      const C& c)
  {
    for(int i=0;i<NDIRAC;i++)
      color_summ_the_prod_double(a[i],b[i],c);
  }
  
  template <typename A,
typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spincolor_prod_idouble(A&& a,
			      const B& b,
			      const double& factor)
  {
    for(int i=0;i<NDIRAC;i++)
      color_prod_idouble(a[i],b[i],factor);
  }
  
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spincolor_prodassign_double(A&& a,
				   const double& factor)
  {
    spincolor_prod_double(a,a,factor);
  }
  
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spincolor_prodassign_idouble(A&& a,
				    const double& factor)
  {
    spincolor_prod_idouble(a,a,factor);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spincolor_prod_color_delta(A&& a,
				  const B& b,
				  const int c)
  {
    for(size_t i=0;i<NDIRAC;i++)
      color_prod_color_delta(a[i],b[i],c);
  }
  
  //spincolor*complex
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_spincolor_prod_complex(A&& a,
				     const B& b,
				     const C& factor)
  {
    for(int i=0;i<NDIRAC;i++)
      unsafe_color_prod_complex(a[i],b[i],factor);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_spincolor_prod_complex(A&& a,
				   const B& b,
				   const C& factor)
  {
    spincolor c;
    spincolor_copy(c,b);
    unsafe_spincolor_prod_complex(a,c,factor);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spincolor_prodassign_complex(A&&a,
				    const B& factor)
  {
    safe_spincolor_prod_complex(a,a,factor);
  }
  
  //spincolor+spincolor*complex
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spincolor_summ_the_prod_complex(A&& a,
				       const B& b,
				       const C& factor)
  {
    for(int i=0;i<NDIRAC;i++)
      color_summ_the_prod_complex(a[i],b[i],factor);
  }
  
  //spincolor+spincolor*idouble
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spincolor_summ_the_prod_idouble(A&& a,
				       const B& b,
				       const C& c)
  {
    for(int i=0;i<NDIRAC;i++)
      color_summ_the_prod_idouble(a[i],b[i],c);
  }
  
  /// Dirac*spincolor
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_dirac_prod_spincolor(A&& out,
				   const dirac_matr& m,
				   const B& in)
  {
    for(int id1=0;id1<NDIRAC;id1++)
      unsafe_color_prod_complex(out[id1],in[m.pos[id1]],m.entr[id1]);
  }
  
  CUDA_HOST_AND_DEVICE inline void dirac_summ_the_prod_spincolor(spincolor out,const dirac_matr& m,const spincolor in)
  {for(size_t id1=0;id1<NDIRAC;id1++) color_summ_the_prod_complex(out[id1],in[m.pos[id1]],m.entr[id1]);}
  
  CUDA_HOST_AND_DEVICE inline void dirac_subt_the_prod_spincolor(spincolor out,const dirac_matr& m,const spincolor in)
  {for(size_t id1=0;id1<NDIRAC;id1++) color_subt_the_prod_complex(out[id1],in[m.pos[id1]],m.entr[id1]);}
  
  //spincolor*dirac
  inline void unsafe_spincolor_prod_dirac(spincolor out,const spincolor in,const dirac_matr& m)
  {spincolor_put_to_zero(out);for(size_t id1=0;id1<NDIRAC;id1++) color_summ_the_prod_complex(out[m.pos[id1]],in[id1],m.entr[id1]);}
  
  //dirac*spincolor
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_dirac_prod_spincolor(A&& out,
				 const dirac_matr& m,
				 const B& in)
  {
    spincolor tmp;
    unsafe_dirac_prod_spincolor(tmp,m,in);
    spincolor_copy(out,tmp);
  }
  
  //spincolor*dirac
  inline void safe_spincolor_prod_dirac(spincolor out,const spincolor in,const dirac_matr& m)
  {spincolor tmp;unsafe_spincolor_prod_dirac(tmp,in,m);spincolor_copy(out,tmp);}
  
  /// su3*spincolor
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3_prod_spincolor(A&& out,
				 const B& U,
				 const C& in)
  {
    for(int is=0;is<NDIRAC;is++)
      unsafe_su3_prod_color(out[is],U,in[is]);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_summ_the_prod_spincolor(A&& out,
				   const B& U,
				   const C& in)
  {
    for(int is=0;is<NDIRAC;is++)
      su3_summ_the_prod_color(out[is],U,in[is]);
  }
  
  CUDA_HOST_AND_DEVICE inline void su3_subt_the_prod_spincolor(spincolor out,const su3 U,const spincolor in)
  {for(size_t is=0;is<NDIRAC;is++) su3_subt_the_prod_color(out[is],U,in[is]);}
  
  /// su3^*spincolor
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3_dag_prod_spincolor(A&& out,
				     const B& U,
				     const C& in)
  {
    for(int is=0;is<NDIRAC;is++)
      unsafe_su3_dag_prod_color(out[is],U,in[is]);
  }
  
  inline void safe_su3_dag_prod_spincolor(spincolor out,const su3 U,const spincolor in)
  {spincolor temp;unsafe_su3_dag_prod_spincolor(temp,U,in);spincolor_copy(out,temp);}
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_dag_summ_the_prod_spincolor(A&& out,
				       const B& U,
				       const C& in)
  {
    for(int is=0;is<NDIRAC;is++)
      su3_dag_summ_the_prod_color(out[is],U,in[is]);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_dag_subt_the_prod_spincolor(A&& out,
				       const B& U,
				       const C& in)
  {
    for(int is=0;is<NDIRAC;is++)
      su3_dag_subt_the_prod_color(out[is],U,in[is]);
  }
  
  /// Product of spincolor and spinspin
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_spincolor_prod_spinspin(A&& a,
				      const B& b,
				      const C& c)
  {
    spincolor_put_to_zero(a);
    
    for(int id1=0;id1<NDIRAC;id1++)
      for(int id2=0;id2<NDIRAC;id2++)
	for(int ic=0;ic<NCOL;ic++)
	  complex_summ_the_prod(a[id1][ic],b[id2][ic],c[id2][id1]);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_spincolor_prod_spinspin(A&& a,
				    const B& b,
				    const C& c)
  {
    spincolor tmp;
    unsafe_spincolor_prod_spinspin(tmp,b,c);
    spincolor_copy(a,tmp);
  }
  
  /// Prouduct of spinspin and spin
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_spinspin_prod_spincolor(A&& a,
				      const B& b,
				      const C& c)
  {
    spincolor_put_to_zero(a);
    for(int id1=0;id1<NDIRAC;id1++)
      for(int id2=0;id2<NDIRAC;id2++)
	for(int ic=0;ic<NCOL;ic++)
	  complex_summ_the_prod(a[id1][ic],b[id1][id2],c[id2][ic]);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_spinspin_prod_spincolor(A&& a,
				    const B& b,
				    const C& c)
  {
    spincolor tmp;
    unsafe_spinspin_prod_spincolor(tmp,b,c);
    spincolor_copy(a,tmp);
  }
  
  //su3^*gamma*spincolor
  inline void unsafe_su3_dag_dirac_prod_spincolor(spincolor out,const su3 U,const dirac_matr& m,const spincolor in)
  {
    color tmp;
    for(size_t id1=0;id1<NDIRAC;id1++)
      {
	for(size_t ic=0;ic<NCOL;ic++) unsafe_complex_prod(tmp[ic],m.entr[id1],in[m.pos[id1]][ic]);
	unsafe_su3_dag_prod_color(out[id1],U,tmp);
      }
  }
  
  inline void unsafe_su3_dag_dirac_summ_the_prod_spincolor(spincolor out,const su3 U,const dirac_matr& m,const spincolor in)
  {
    color tmp;
    for(size_t id1=0;id1<NDIRAC;id1++)
      {
	for(size_t ic=0;ic<NCOL;ic++) unsafe_complex_prod(tmp[ic],m.entr[id1],in[m.pos[id1]][ic]);
	su3_dag_summ_the_prod_color(out[id1],U,tmp);
      }
  }
  
  //su3*dirac*spincolor
  inline void unsafe_su3_dirac_prod_spincolor(spincolor out,const su3 U,const dirac_matr& m,const spincolor in)
  {
    color tmp;
    for(size_t id1=0;id1<NDIRAC;id1++)
      {
	for(size_t ic=0;ic<NCOL;ic++) unsafe_complex_prod(tmp[ic],m.entr[id1],in[m.pos[id1]][ic]);
	unsafe_su3_prod_color(out[id1],U,tmp);
      }
  }
  
  inline void unsafe_su3_dirac_subt_the_prod_spincolor(spincolor out,const su3 U,const dirac_matr& m,const spincolor in)
  {
    color tmp;
    for(size_t id1=0;id1<NDIRAC;id1++)
      {
	for(size_t ic=0;ic<NCOL;ic++) unsafe_complex_prod(tmp[ic],m.entr[id1],in[m.pos[id1]][ic]);
	su3_subt_the_prod_color(out[id1],U,tmp);
      }
  }
  
  inline void unsafe_su3_dirac_summ_the_prod_spincolor(spincolor out,const su3 U,const dirac_matr& m,const spincolor in)
  {
    color tmp;
    for(size_t id1=0;id1<NDIRAC;id1++)
      {
	for(size_t ic=0;ic<NCOL;ic++) unsafe_complex_prod(tmp[ic],m.entr[id1],in[m.pos[id1]][ic]);
	su3_summ_the_prod_color(out[id1],U,tmp);
      }
  }
  
  /// Take the complex conjugated of the first
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spincolor_scalar_prod(A&& out,
			     const B& in1,
			     const C& in2)
  {
    complex_put_to_zero(out);
    for(int id=0;id<NDIRAC;id++)
      {
	complex temp;
	color_scalar_prod(temp,in1[id],in2[id]);
	complex_summassign(out,temp);
      }
  }
  
  inline double spincolor_norm2(const spincolor s)
  {double out=color_norm2(s[0]);for(size_t id=1;id<NDIRAC;id++) out+=color_norm2(s[id]);return out;}
  
  ///////////////////////////////// su3*colorspinspin ///////////////////////////////////
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3_prod_colorspinspin(A&& a,
				     const B& b,
				     const C& c)
  {
    for(int id_so=0;id_so<NDIRAC;id_so++)
      for(int id_si=0;id_si<NDIRAC;id_si++)
	for(int c1=0;c1<NCOL;c1++)
	  {
	    unsafe_complex_prod(a[c1][id_si][id_so],b[c1][0],c[0][id_si][id_so]);
	    for(int c2=1;c2<NCOL;c2++) complex_summ_the_prod(a[c1][id_si][id_so],b[c1][c2],c[c2][id_si][id_so]);
	  }
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3_dag_prod_colorspinspin(A&& a,
					 const B& b,
					 const C& c)
  {
    for(int id_so=0;id_so<NDIRAC;id_so++)
      for(int id_si=0;id_si<NDIRAC;id_si++)
	for(int c1=0;c1<NCOL;c1++)
	  {
	    unsafe_complex_conj1_prod(a[c1][id_si][id_so],b[0][c1],c[0][id_si][id_so]);
	    for(int c2=1;c2<NCOL;c2++)
	      complex_summ_the_conj1_prod(a[c1][id_si][id_so],b[c2][c1],c[c2][id_si][id_so]);
	  }
  }
  
  CUDA_HOST_AND_DEVICE inline void su3_summ_the_prod_colorspinspin(colorspinspin a,const su3 b,const colorspinspin c)
  {
    for(size_t id_so=0;id_so<NDIRAC;id_so++)
      for(size_t id_si=0;id_si<NDIRAC;id_si++)
	for(size_t c1=0;c1<NCOL;c1++)
	  for(size_t c2=0;c2<NCOL;c2++)
	    complex_summ_the_prod(a[c1][id_si][id_so],b[c1][c2],c[c2][id_si][id_so]);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_dag_subt_the_prod_colorspinspin(A&& a,
					   const B& b,
					   const C& c)
  {
    for(int id_so=0;id_so<NDIRAC;id_so++)
      for(int id_si=0;id_si<NDIRAC;id_si++)
	for(int c1=0;c1<NCOL;c1++)
	  for(int c2=0;c2<NCOL;c2++)
	    complex_subt_the_conj1_prod(a[c1][id_si][id_so],b[c2][c1],c[c2][id_si][id_so]);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_dag_summ_the_prod_colorspinspin(A&& a,
					   const B& b,
					   const C& c)
  {
    for(int id_so=0;id_so<NDIRAC;id_so++)
      for(int id_si=0;id_si<NDIRAC;id_si++)
	for(int c1=0;c1<NCOL;c1++)
	  for(int c2=0;c2<NCOL;c2++)
	    complex_summ_the_conj1_prod(a[c1][id_si][id_so],b[c2][c1],c[c2][id_si][id_so]);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_dirac_prod_colorspinspin(A&&out,
				       const dirac_matr& m,
				       const B& in)
  {
    for(int ic=0;ic<NCOL;ic++)
      unsafe_dirac_prod_spinspin(out[ic],m,in[ic]);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_dirac_prod_colorspinspin(A&&out,
				     const dirac_matr& m,
				     const B& in)
  {
    colorspinspin temp;
    unsafe_dirac_prod_colorspinspin(temp,m,in);
    colorspinspin_copy(out,temp);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_dirac_prod_su3spinspin(A&&out,
				     const dirac_matr& m,
				     const B& in)
  {
    for(int ic=0;ic<NCOL;ic++)
      unsafe_dirac_prod_colorspinspin(out[ic],m,in[ic]);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_dirac_prod_su3spinspin(A&& out,
				   const dirac_matr& m,
				   const B& in)
  {
    su3spinspin temp;
    unsafe_dirac_prod_su3spinspin(temp,m,in);
    su3spinspin_copy(out,temp);
  }
  
  //////////////////////////////// get and put ///////////////////////////////////
  
  //Get a color from a colorspinspin
  CUDA_HOST_AND_DEVICE inline void get_color_from_colorspinspin(color out,const colorspinspin in,int id1,int id2)
  {for(int ic_sink=0;ic_sink<NCOL;ic_sink++) complex_copy(out[ic_sink],in[ic_sink][id1][id2]);}
  
  //Get a color from a spincolor
  CUDA_HOST_AND_DEVICE inline void get_color_from_spincolor(color out,const spincolor in,int id)
  {for(int ic_sink=0;ic_sink<NCOL;ic_sink++) complex_copy(out[ic_sink],in[id][ic_sink]);}
  
  //Get a spincolor from a colorspinspin
  //In a spinspin the sink index runs slower than the source
  CUDA_HOST_AND_DEVICE inline void get_spincolor_from_colorspinspin(spincolor out,const colorspinspin in,int id_source)
  {
    for(int ic_sink=0;ic_sink<NCOL;ic_sink++)
      for(int id_sink=0;id_sink<NDIRAC;id_sink++) //dirac index of sink
	complex_copy(out[id_sink][ic_sink],in[ic_sink][id_sink][id_source]);
  }
  
  //Get a spincolor from a su3spinspin
  CUDA_HOST_AND_DEVICE inline void get_spincolor_from_su3spinspin(spincolor out,const su3spinspin in,int id_source,int ic_source)
  {
    for(int ic_sink=0;ic_sink<NCOL;ic_sink++)
      for(int id_sink=0;id_sink<NDIRAC;id_sink++) //dirac index of sink
	complex_copy(out[id_sink][ic_sink],in[ic_sink][ic_source][id_sink][id_source]);
  }
  
  //Get a color from a su3
  CUDA_HOST_AND_DEVICE inline void get_color_from_su3(color out,const su3 in,int ic_source)
  {for(int ic_sink=0;ic_sink<NCOL;ic_sink++) complex_copy(out[ic_sink],in[ic_sink][ic_source]);}
  
  //Put a color into a colorspinspin
  CUDA_HOST_AND_DEVICE inline void put_color_into_colorspinspin(colorspinspin out,color in,int id1,int id2)
  {for(int ic_sink=0;ic_sink<NCOL;ic_sink++) complex_copy(out[ic_sink][id1][id2],in[ic_sink]);}
  
  //Put a color into a spincolor
  CUDA_HOST_AND_DEVICE inline void put_color_into_spincolor(spincolor out,color in,int id)
  {for(int ic_sink=0;ic_sink<NCOL;ic_sink++) complex_copy(out[id][ic_sink],in[ic_sink]);}
  
  //Put a spincolor into a colorspinspin
  CUDA_HOST_AND_DEVICE inline void put_spincolor_into_colorspinspin(colorspinspin out,const spincolor in,int id_source)
  {
    for(int ic_sink=0;ic_sink<NCOL;ic_sink++)
      for(int id_sink=0;id_sink<NDIRAC;id_sink++) //dirac index of sink
	complex_copy(out[ic_sink][id_sink][id_source],in[id_sink][ic_sink]);
  }
  
  //Put a spincolor into a su3spinspin
  CUDA_HOST_AND_DEVICE inline void put_spincolor_into_su3spinspin(su3spinspin out,const spincolor in,int id_source,int ic_source)
  {
    for(int ic_sink=0;ic_sink<NCOL;ic_sink++)
      for(int id_sink=0;id_sink<NDIRAC;id_sink++) //dirac index of sink
	complex_copy(out[ic_sink][ic_source][id_sink][id_source],in[id_sink][ic_sink]);
  }
  
  //Put a spincolor into a su3
  CUDA_HOST_AND_DEVICE inline void put_color_into_su3(su3 out,color in,int ic_source)
  {for(int ic_sink=0;ic_sink<NCOL;ic_sink++) complex_copy(out[ic_sink][ic_source],in[ic_sink]);}
  ///////////////////////////////////// colorspinspin ////////////////////////////////////
  
  //colorspinspin*real or complex
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void colorspinspin_prod_double(A&& out,
				 const B& in,
				 const double& factor)
  {
    for(size_t i=0;i<NCOL;i++)
      spinspin_prod_double(out[i],in[i],factor);
  }
  
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void colorspinspin_prodassign_double(A&& c,
				       const double& f)
  {
    colorspinspin_prod_double(c,c,f);
  }
  
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void colorspinspin_prodassign_idouble(A&& c,
					const double& f)
  {
    for(int ic=0;ic<NCOL;ic++)
      spinspin_prodassign_idouble(c[ic],f);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void colorspinspin_prod_complex(A&& out,
				  const B& in,
				  const C& factor)
  {
    for(size_t ic=0;ic<NCOL;ic++)
      safe_spinspin_prod_complex(out[ic],in[ic],factor);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void colorspinspin_prod_complex_conj(A&& out,
								   const B& in,
								   const C& factor)
  {
    complex temp;
    complex_conj(temp,factor);
    colorspinspin_prod_complex(out,in,temp);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void colorspinspin_prodassign_complex(A&& c,
					const B& f)
  {
    colorspinspin_prod_complex(c,c,f);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void colorspinspin_prodassign_complex_conj(A&& c,
					     const B& f)
  {
    colorspinspin_prod_complex_conj(c,c,f);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void colorspinspin_prod_color_delta(A&& a,
				      const B& b,
				      const size_t& c)
  {
    for(size_t ic=0;ic<NCOL;ic++)
      spinspin_prod_double(a[ic],b[ic],c==ic);
  }
  
  //colorspinspin summ
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void colorspinspin_summ(A&& out,
			  const B& in1,
			  const C& in2)
  {
    for(int i=0;i<NCOL;i++)
      spinspin_summ(out[i],in1[i],in2[i]);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void colorspinspin_summassign(A&& out,
				const B& in)
  {
    colorspinspin_summ(out,out,in);
  }
  
  //colorspinspin subt
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void colorspinspin_subt(A&& out,
			  const B& in1,
			  const C& in2)
  {
    for(size_t i=0;i<NCOL;i++)
      spinspin_subt(out[i],in1[i],in2[i]);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void colorspinspin_subtassign(A&& out,
				const B& in)
  {
    colorspinspin_subt(out,out,in);
  }
  
  //summ two colorspinspin with a factor
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void colorspinspin_summ_the_prod_double(A&& out,
					  const B& in,
					  const double& factor)
  {
    for(int i=0;i<NCOL;i++)
      spinspin_summ_the_prod_double(out[i],in[i],factor);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void colorspinspin_summ_the_prod_idouble(A&& out,
					   const B& in,
					   const double& factor)
  {
    for(size_t i=0;i<NCOL;i++)
      spinspin_summ_the_prod_idouble(out[i],in[i],factor);
  }
  
  //colorspinspin*complex
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_colorspinspin_prod_complex(A&& out,
					 const B& in,
					 const C& factor)
  {
    for(size_t i=0;i<NCOL;i++)
      unsafe_spinspin_prod_complex(out[i],in[i],factor);
  }
  
  /////////////////////////////// su3spinspin /////////////////////////////////////////////
  
  //colorspinspin*real or complex
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3spinspin_prod_double(A&& out,
			       const B& in,
			       const double& factor)
  {
    for(size_t i=0;i<NCOL;i++)
      colorspinspin_prod_double(out[i],in[i],factor);
  }
  
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3spinspin_prodassign_double(A&& out,
				     const double& factor)
  {
    for(int i=0;i<NCOL;i++)
      colorspinspin_prodassign_double(out[i],factor);
  }
  
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3spinspin_prodassign_idouble(A&& out,
				      const double& factor)
  {
    for(int i=0;i<NCOL;i++)
      colorspinspin_prodassign_idouble(out[i],factor);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3spinspin_prod_complex(A&& a,
				const B& b,
				const C& c)
  {
    for(size_t i=0;i<NCOL;i++)
      colorspinspin_prod_complex(a[i],b[i],c);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3spinspin_prod_complex_conj(A&& a,
				     const B& b,
				     const C& c)
  {
    for(size_t i=0;i<NCOL;i++)
      colorspinspin_prod_complex_conj(a[i],b[i],c);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3spinspin_prodassign_complex(A&& a,
				      const B& b)
  {
    for(size_t i=0;i<NCOL;i++)
      colorspinspin_prodassign_complex(a[i],b);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3spinspin_prodassign_complex_conj(A&& a,
					   const B& b)
  {
    for(size_t i=0;i<NCOL;i++)
      colorspinspin_prodassign_complex_conj(a[i],b);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3spinspin_prod_color_delta(A&& a,
				    const B& b,
				    const size_t& c)
  {
    for(size_t ic=0;ic<NCOL;ic++)
      colorspinspin_prod_double(a[ic],b[ic],c==ic);
  }
  
  //su3spinspin summ
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3spinspin_summ(A&& out,
			const B& in1,
			const C& in2)
  {
    for(int i=0;i<NCOL;i++)
      colorspinspin_summ(out[i],in1[i],in2[i]);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3spinspin_summassign(A&& out,
			      const B& in)
  {
    su3spinspin_summ(out,out,in);
  }
  
  //su3spinspin subt
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3spinspin_subt(A&& a,
			const B& b,
			const C& c)
  {
    for(size_t i=0;i<NCOL;i++)
      colorspinspin_subt(a[i],b[i],c[i]);}
  
  //su3spinspin subt
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3spinspin_subtassign(A&& a,
			      const B& b)
  {
    su3spinspin_subt(a,a,b);
  }
  
  //summ two su3spinspin with a factor
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3spinspin_summ_the_prod_double(A&& out,
					const B& in,
					const double& factor)
  {
    for(int i=0;i<NCOL;i++)
      colorspinspin_summ_the_prod_double(out[i],in[i],factor);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3spinspin_summ_the_prod_idouble(A&& out,
					 const B& in,
					 const double& factor)
  {
    for(size_t i=0;i<NCOL;i++)
      colorspinspin_summ_the_prod_idouble(out[i],in[i],factor);
  }
  
  ////////////////////////////////////// su3spinspin /////////////////////////////////////
  
  //su3spinspin*complex
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3spinspin_prod_complex(A&& out,
				       const B& in,
				       const C& factor)
  {
    for(size_t i=0;i<NCOL;i++)
      unsafe_colorspinspin_prod_complex(out[i],in[i],factor);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_su3spinspin_prod_complex(A&& out,
				     const B& in,
				     const C& factor)
  {
    su3spinspin temp;
    su3spinspin_copy(temp,in);
    unsafe_su3spinspin_prod_complex(out,temp,factor);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3_prod_su3spinspin(A&& a,
				   const B& b,
				   const C& c)
  {
    for(int id_so=0;id_so<NDIRAC;id_so++)
      for(int id_si=0;id_si<NDIRAC;id_si++)
	for(int c1=0;c1<NCOL;c1++)
	  for(int c3=0;c3<NCOL;c3++)
	    {
	      unsafe_complex_prod(a[c1][c3][id_si][id_so],b[c1][0],c[0][c3][id_si][id_so]);
	      for(int c2=1;c2<NCOL;c2++)
		complex_summ_the_prod(a[c1][c3][id_si][id_so],b[c1][c2],c[c2][c3][id_si][id_so]);
	    }
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_su3_dag_prod_su3spinspin(A&& a,
				       const B& b,
				       const C& c)
  {
    for(int id_so=0;id_so<NDIRAC;id_so++)
      for(int id_si=0;id_si<NDIRAC;id_si++)
	for(int c1=0;c1<NCOL;c1++)
	  for(int c3=0;c3<NCOL;c3++)
	    {
	      unsafe_complex_conj1_prod(a[c1][c3][id_si][id_so],b[0][c1],c[0][c3][id_si][id_so]);
	      for(int c2=1;c2<NCOL;c2++)
		complex_summ_the_conj1_prod(a[c1][c3][id_si][id_so],b[c2][c1],c[c2][c3][id_si][id_so]);
	    }
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_summ_the_prod_su3spinspin(A&& a,
				     const B& b,
				     const C& c)
  {
    for(size_t c1=0;c1<NCOL;c1++)
      for(size_t c3=0;c3<NCOL;c3++)
	for(size_t c2=0;c2<NCOL;c2++)
	  for(size_t id_si=0;id_si<NDIRAC;id_si++)
	    for(size_t id_so=0;id_so<NDIRAC;id_so++)
	      complex_summ_the_prod(a[c1][c3][id_si][id_so],b[c1][c2],c[c2][c3][id_si][id_so]);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_dag_subt_the_prod_su3spinspin(A&& a,
					 const B& b,
					 const C& c)
  {
    for(int c1=0;c1<NCOL;c1++)
      for(int c3=0;c3<NCOL;c3++)
	for(int c2=0;c2<NCOL;c2++)
	  for(int id_si=0;id_si<NDIRAC;id_si++)
	    for(int id_so=0;id_so<NDIRAC;id_so++)
	      complex_subt_the_conj1_prod(a[c1][c3][id_si][id_so],b[c2][c1],c[c2][c3][id_si][id_so]);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void su3_dag_summ_the_prod_su3spinspin(A&& a,
					 const B& b,
					 const C& c)
  {
    for(int id_so=0;id_so<NDIRAC;id_so++)
      for(int id_si=0;id_si<NDIRAC;id_si++)
	for(int c1=0;c1<NCOL;c1++)
	  for(int c2=0;c2<NCOL;c2++)
	    for(int c3=0;c3<NCOL;c3++)
	      complex_summ_the_conj1_prod(a[c1][c3][id_si][id_so],b[c2][c1],c[c2][c3][id_si][id_so]);
  }
  
  CUDA_HOST_AND_DEVICE void hermitian_exact_i_exponentiate_ingredients(hermitian_exp_ingredients &out,
								       const su3& Q);
  
  /// Build the exponential from the ingredients
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_hermitian_exact_i_exponentiate(A&& out,
					   const hermitian_exp_ingredients& ing)
  {
    //compute out according to (eq. 13)
    su3_put_to_diag_complex(out,ing.f[0]);
    su3_summ_the_prod_complex(out,ing.Q,ing.f[1]);
    su3_summ_the_prod_complex(out,ing.Q2,ing.f[2]);
  }
  
  template <typename A,
	    ENABLE_THIS_TEMPLATE_IF(not std::is_same_v<A,su3>)>
  CUDA_HOST_AND_DEVICE void hermitian_exact_i_exponentiate_ingredients(hermitian_exp_ingredients &out,
								       const A& Q)
  {
    su3 temp;
    su3_copy(temp,Q);
    hermitian_exact_i_exponentiate_ingredients(out,temp);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE inline void safe_hermitian_exact_i_exponentiate(A&& out,
								       const B& Q)
  {
    hermitian_exp_ingredients ing;
    hermitian_exact_i_exponentiate_ingredients(ing,Q);
    safe_hermitian_exact_i_exponentiate(out,ing);
  }
  
  template <typename A,
	    typename B>
  /// Can be used directly from an anti-hermitian matrix, ie. it compute straight exp(iQ)
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_anti_hermitian_exact_exponentiate(A&& out,
					      const B& iQ)
  {
    su3 Q;
    su3_prod_idouble(Q,iQ,-1);
    
    safe_hermitian_exact_i_exponentiate(out,Q);
  }
  
  //return sqrt(|U*U^+-1|)
  template <typename U>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  double su3_get_non_unitariness(const U& u)
  {
    su3 zero;
    su3_put_to_id(zero);
    su3_subt_the_prod_su3_dag(zero,u,u);
    
    return
      sqrt(real_part_of_trace_su3_prod_su3_dag(zero,zero));
  }
}

#endif
