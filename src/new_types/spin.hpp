#ifndef _SPIN_HPP
#define _SPIN_HPP

#include <string.h>
#include <stdio.h>

#include "complex.hpp"
#include "dirac.hpp"

#include "routines/math_routines.hpp"

#ifndef EXTERN_SPIN
 #define EXTERN_SPIN extern
#endif

namespace nissa
{
  typedef complex spin[NDIRAC];
  typedef complex halfspin[2];
  
  typedef spin spinspin[NDIRAC];
  typedef complex as2t[NDIM*(NDIM+1)/2];
  
  //this is just for avoid misleading, but is nothing more that a spinspin
  typedef complex spin1field[NDIM];
  typedef spin1field spin1prop[NDIM];
  
  EXTERN_SPIN as2t smunu_entr[4];   //these are the sigma matrices entries
  EXTERN_SPIN int smunu_pos[4][6];  //and positions
  EXTERN_SPIN spinspin opg[4],omg[4];
  
  /////////////////////////////////////////////////////
  
  //Print a spin
  inline void spin_print(const spin s)
  {
    for(int id=0;id<NDIRAC;id++) printf("%+16.16lg,%+16.16lg\t",s[id][0],s[id][1]);
    printf("\n");
  }
  
  /// s=0
  template <typename S>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spin_put_to_zero(S&& s)
  {
    for(int i=0;i<NDIRAC;i++)
      complex_put_to_zero(s[i]);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spin_copy(A&& out,
		 const B& in)
  {
    for(int i=0;i<NDIRAC;i++)
      complex_copy(out[i],in[i]);
  }
  
  inline void spin_conj(spin out,const spin in)
  {for(int i=0;i<NDIRAC;i++) complex_conj(out[i],in[i]);}
  
  inline void spin_summ(spin a,const spin b,const spin c)
  {for(int i=0;i<NDIRAC;i++) complex_summ(a[i],b[i],c[i]);}
  inline void spin_summassign(spin a,const spin b)
  {spin_summ(a,a,b);}
  inline void spin_subt(spin a,const spin b,const spin c)
  {for(int i=0;i<NDIRAC;i++) complex_subt(a[i],b[i],c[i]);}
  inline void spin_subtassign(spin a,const spin b)
  {spin_subt(a,a,b);}
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spin_prod_double(A&& a,
			const B& b,
			const double& c)
  {
    for(int i=0;i<NDIRAC;i++)
      complex_prod_double(a[i],b[i],c);
  }
  
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spin_prodassign_double(A&& a,
			      const double& b)
  {
    spin_prod_double(a,a,b);
  }
  
  inline void spin_summ_the_complex_prod(spin a,const spin b,const complex c)
  {for(int i=0;i<NDIRAC;i++) complex_summ_the_prod(a[i],b[i],c);}
  inline void spin_subt_the_complex_prod(spin a,const spin b,const complex c)
  {for(int i=0;i<NDIRAC;i++) complex_subt_the_prod(a[i],b[i],c);}
  
  inline void spin_summ_the_complex_conj2_prod(spin a,const spin b,const complex c)
  {for(int i=0;i<NDIRAC;i++) complex_summ_the_conj2_prod(a[i],b[i],c);}
  inline void spin_subt_the_complex_conj2_prod(spin a,const spin b,const complex c)
  {for(int i=0;i<NDIRAC;i++) complex_subt_the_conj2_prod(a[i],b[i],c);}
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spinspin_copy(A&& a,
		     const B& b)
  {
    for(int i=0;i<NDIRAC;i++)
      spin_copy(a[i],b[i]);
  }
  
  inline void unsafe_spinspin_hermitian(spinspin b,const spinspin a)
  {
    for(int id1=0;id1<NDIRAC;id1++)
      for(int id2=0;id2<NDIRAC;id2++)
	complex_conj(b[id1][id2],a[id2][id1]);
  }
  inline void safe_spinspin_hermitian(spinspin b,const spinspin a)
  {
    spinspin temp;
    unsafe_spinspin_hermitian(temp,a);
    spinspin_copy(b,temp);
  }
  
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spinspin_put_to_zero(A&& a)
  {
    for(int i=0;i<NDIRAC;i++)
      spin_put_to_zero(a[i]);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spinspin_put_to_diag_complex(A&& a,
				    const B& b)
  {
    spinspin_put_to_zero(a);
    for(int id=0;id<NDIRAC;id++)
      complex_copy(a[id][id],b);
  }
  
  template <typename A>
  inline void spinspin_put_to_diag_double(A&& a,
					  const double& b)
  {
    const complex c={b,0};
    
    spinspin_put_to_diag_complex(a,c);
  }
  
  inline void spinspin_put_to_idiag_double(spinspin a,double b)
  {const complex c={0,b};spinspin_put_to_diag_complex(a,c);}
  
  inline void spin_direct_prod(spinspin out,const spin a,const spin b)
  {
    for(int id=0;id<NDIRAC;id++)
      for(int jd=0;jd<NDIRAC;jd++)
	unsafe_complex_prod(out[id][jd],a[id],b[jd]);
  }
  
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spinspin_put_to_id(A&& a)
  {
    const complex o={1,0};
    
    spinspin_put_to_diag_complex(a,o);
  }
  
  CUDA_HOST_AND_DEVICE inline void spinspin_summ(spinspin a,const spinspin b,const spinspin c)
  {for(int id1=0;id1<NDIRAC;id1++) for(int id2=0;id2<NDIRAC;id2++) complex_summ(a[id1][id2],b[id1][id2],c[id1][id2]);}
  CUDA_HOST_AND_DEVICE inline void spinspin_subt(spinspin a,const spinspin b,const spinspin c)
  {for(int id1=0;id1<NDIRAC;id1++) for(int id2=0;id2<NDIRAC;id2++) complex_subt(a[id1][id2],b[id1][id2],c[id1][id2]);}
  CUDA_HOST_AND_DEVICE inline void spinspin_summassign(spinspin a,const spinspin b)
  {spinspin_summ(a,a,b);}
  inline void spinspin_subtassign(spinspin a,const spinspin b)
  {spinspin_subt(a,a,b);}
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE inline void spinspin_prod_double(A&& a,
							const B& b,
							const double& c)
  {
    for(int id1=0;id1<NDIRAC;id1++)
      for(int id2=0;id2<NDIRAC;id2++)
	complex_prod_double(a[id1][id2],b[id1][id2],c);
  }
  
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spinspin_prodassign_double(A&& a,
				  const double& b)
  {
    spinspin_prod_double(a,a,b);
  }
  
  CUDA_HOST_AND_DEVICE inline void spinspin_summ_the_prod_double(spinspin a,const spinspin b,double c)
  {for(int id1=0;id1<NDIRAC;id1++) for(int id2=0;id2<NDIRAC;id2++) complex_summ_the_prod_double(a[id1][id2],b[id1][id2],c);}
  CUDA_HOST_AND_DEVICE inline void spinspin_prod_idouble(spinspin a,const spinspin b,double c)
  {for(int id1=0;id1<NDIRAC;id1++) for(int id2=0;id2<NDIRAC;id2++) complex_prod_idouble(a[id1][id2],b[id1][id2],c);}
  CUDA_HOST_AND_DEVICE inline void spinspin_prodassign_idouble(spinspin a,double b)
  {spinspin_prod_idouble(a,a,b);}
  CUDA_HOST_AND_DEVICE inline void spinspin_summ_the_prod_idouble(spinspin a,const spinspin b,double c)
  {for(int id1=0;id1<NDIRAC;id1++) for(int id2=0;id2<NDIRAC;id2++) complex_summ_the_prod_idouble(a[id1][id2],b[id1][id2],c);}
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spinspin_summ_the_complex_prod(A&& a,
				      const B& b,
				      const C& c)
  {
    for(int id1=0;id1<NDIRAC;id1++)
      for(int id2=0;id2<NDIRAC;id2++)
	complex_summ_the_prod(a[id1][id2],b[id1][id2],c);
  }
  
  inline void spinspin_subt_the_complex_prod(spinspin a,const spinspin b,const complex c)
  {for(int id1=0;id1<NDIRAC;id1++) for(int id2=0;id2<NDIRAC;id2++) complex_subt_the_prod(a[id1][id2],b[id1][id2],c);}
  inline void spinspin_summassign_the_complex_prod(spinspin a,const complex b){spinspin_summ_the_complex_prod(a,a,b);}
  inline void spinspin_subtassign_the_complex_prod(spinspin a,const complex b){spinspin_subt_the_complex_prod(a,a,b);}
  
  inline void spinspin_summ_the_complex_conj2_prod(spinspin a,const spinspin b,const complex c)
  {for(int id1=0;id1<NDIRAC;id1++) for(int id2=0;id2<NDIRAC;id2++) complex_summ_the_conj2_prod(a[id1][id2],b[id1][id2],c);}
  inline void spinspin_subt_the_complex_conj2_prod(spinspin a,const spinspin b,const complex c)
  {for(int id1=0;id1<NDIRAC;id1++) for(int id2=0;id2<NDIRAC;id2++) complex_subt_the_conj2_prod(a[id1][id2],b[id1][id2],c);}
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_spinspin_prod_complex(A&& a,
				    const B& b,
				    const C& c)
  {
    spinspin_put_to_zero(a);
    spinspin_summ_the_complex_prod(a,b,c);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_spinspin_prod_complex(A&& a,
				  const B& b,
				  const C& c)
  {
    spinspin d;
    unsafe_spinspin_prod_complex(d,b,c);
    spinspin_copy(a,d);
  }
  
  CUDA_HOST_AND_DEVICE inline void unsafe_spinspin_prod_complex_conj2(spinspin a,const spinspin b,const complex c)
  {
    complex d;
    complex_conj(d,c);
    unsafe_spinspin_prod_complex(a,b,d);
  }
  CUDA_HOST_AND_DEVICE inline void safe_spinspin_prod_complex_conj2(spinspin a,const spinspin b,const complex c)
  {
    complex d;
    complex_conj(d,c);
    safe_spinspin_prod_complex(a,b,d);
  }
  
  //saturate two anti-simmetric tensors
  inline void as2t_saturate(complex out,as2t a,as2t b)
  {
    unsafe_complex_prod(out,a[0],b[0]);
    for(int munu=1;munu<6;munu++) complex_summ_the_prod(out,a[munu],b[munu]);
  }
  
  //Print a spinspin
  inline void spinspin_print(const spinspin s)
  {
    for(int id1=0;id1<NDIRAC;id1++)
      {
	for(int id2=0;id2<NDIRAC;id2++) printf("%+16.16lg,%+16.16lg\t",s[id1][id2][0],s[id1][id2][1]);
	printf("\n");
      }
  }
  
  //trace of the product with a dirac matr of a spinspin
  inline void summ_the_trace_dirac_prod_spinspin(complex c,const dirac_matr& a,const spinspin b)
  {
    for(int id=0;id<NDIRAC;id++)
      complex_summ_the_prod(c,a.entr[id],b[a.pos[id]][id]);
  }
  inline void trace_dirac_prod_spinspin(complex c,const dirac_matr& a,const spinspin b)
  {
    c[0]=c[1]=0;
    summ_the_trace_dirac_prod_spinspin(c,a,b);
  }
  
  inline void summ_the_trace_spinspin(complex c,const spinspin a)
  {for(int id=0;id<NDIRAC;id++) complex_summassign(c,a[id][id]);}
  inline void trace_spinspin(complex c,const spinspin a)
  {
    c[0]=c[1]=0;
    summ_the_trace_spinspin(c,a);
  }
  
  //product of two spinspins
  CUDA_HOST_AND_DEVICE inline void spinspin_summ_the_spinspin_dag_prod(spinspin out,const spinspin a,const spinspin b)
  {
    for(int id1=0;id1<NDIRAC;id1++)
      for(int id2=0;id2<NDIRAC;id2++)
	for(int id=0;id<NDIRAC;id++)
	  complex_summ_the_conj2_prod(out[id1][id2],a[id1][id],b[id2][id]);
  }
  CUDA_HOST_AND_DEVICE inline void unsafe_spinspin_prod_spinspin_dag(spinspin out,const spinspin a,const spinspin b)
  {
    memset(out,0,sizeof(spinspin));
    spinspin_summ_the_spinspin_dag_prod(out,a,b);
  }
  inline void safe_spinspin_prod_spinspin_dag(spinspin out,const spinspin a,const spinspin b)
  {
    spinspin c;
    unsafe_spinspin_prod_spinspin_dag(c,a,b);
    memcpy(out,c,sizeof(spinspin));
  }

  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spinspin_summ_the_spinspin_prod(A&& out,
				       const B& a,
				       const C& b)
  {
    for(int id1=0;id1<NDIRAC;id1++)
      for(int id2=0;id2<NDIRAC;id2++)
	for(int id=0;id<NDIRAC;id++)
	  complex_summ_the_prod(out[id1][id2],a[id1][id],b[id][id2]);
  }
  
  inline void spinspin_subt_the_spinspin_prod(spinspin out,const spinspin a,const spinspin b)
  {
    for(int id1=0;id1<NDIRAC;id1++)
      for(int id2=0;id2<NDIRAC;id2++)
	for(int id=0;id<NDIRAC;id++)
	  complex_subt_the_prod(out[id1][id2],a[id1][id],b[id][id2]);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE inline void unsafe_spinspin_prod_spinspin(A&& out,
								 const B& a,
								 const C& b)
  {
    spinspin_put_to_zero(out);
    spinspin_summ_the_spinspin_prod(out,a,b);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_spinspin_prod_spinspin(A&& out,
				   const B& a,
				   const C& b)
  {
    spinspin c;
    unsafe_spinspin_prod_spinspin(c,a,b);
    spinspin_copy(out,c);
  }
  
  inline double real_part_of_trace_spinspin_prod_spinspin_dag(const spinspin a,const spinspin b)
  {
    double t=0;
    
    for(int id1=0;id1<NDIRAC;id1++)
      for(int id2=0;id2<NDIRAC;id2++)
	t+=a[id1][id2][0]*b[id1][id2][0]+a[id1][id2][1]*b[id1][id2][1];
    
    return t;
  }
  
  inline double spinspin_norm2(const spinspin a)
  {return real_part_of_trace_spinspin_prod_spinspin_dag(a,a);}
  
  inline void trace_spinspin_prod_spinspin(complex c,const spinspin a,const spinspin b)
  {
    complex_put_to_zero(c);
    for(int id1=0;id1<NDIRAC;id1++)
      for(int id2=0;id2<NDIRAC;id2++)
	complex_summ_the_prod(c,a[id1][id2],b[id2][id1]);
  }
  
  //Summ the passed gamma multiplied by a double to spinspin
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spinspin_dirac_summ_the_prod_double(A&& out,
					   const dirac_matr& in,
					   const double& r)
  {
    //This is the line on the matrix
    for(int ig=0;ig<NDIRAC;ig++)
      complex_summ_the_prod_double(out[ig][in.pos[ig]],in.entr[ig],r);
  }
  
  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spinspin_dirac_summ_the_prod_idouble(A&& out,
					    const dirac_matr& in,
					    const double& r)
  {
    for(int ig=0;ig<NDIRAC;ig++)
      complex_summ_the_prod_idouble(out[ig][in.pos[ig]],in.entr[ig],r);
  }
  
  template <typename A,
	    typename B>
  CUDA_DEVICE INLINE_FUNCTION
  void spinspin_dirac_summ_the_prod_complex(A&& out,
					    const dirac_matr& in,
					    const B& c)
  {
    for(int ig=0;ig<NDIRAC;ig++)
      complex_summ_the_prod(out[ig][in.pos[ig]],in.entr[ig],c);
  }
  
  inline void spinspin_dirac_subt_the_prod_complex(spinspin out,const dirac_matr& in,const complex c)
  {
    for(int ig=0;ig<NDIRAC;ig++)
      complex_subt_the_prod(out[ig][in.pos[ig]],in.entr[ig],c);
  }

  template <typename A>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spinspin_dirac_prod_double(A&& out,
				  const dirac_matr& in,
				  const double& r)
  {
    spinspin_put_to_zero(out);
    spinspin_dirac_summ_the_prod_double(out,in,r);
  }
  
  inline void spinspin_dirac_prod_idouble(spinspin out,const dirac_matr&in,double r)
  {spinspin_put_to_zero(out);spinspin_dirac_summ_the_prod_idouble(out,in,r);}
  inline void spinspin_dirac_prod_complex(spinspin out,const dirac_matr& in,const complex c)
  {spinspin_put_to_zero(out);spinspin_dirac_summ_the_prod_complex(out,in,c);}
  
  //out=m*in
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_dirac_prod_spinspin(A&& out,
				  const dirac_matr& m,
				  const B& in)
  {
    for(int id1=0;id1<NDIRAC;id1++)
      for(int id2=0;id2<NDIRAC;id2++)
	unsafe_complex_prod(out[id1][id2],m.entr[id1],in[m.pos[id1]][id2]);
  }
  
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_dirac_prod_spinspin(A&& out,
				  const dirac_matr& m,
				  const B& in)
  {
    spinspin temp;
    unsafe_dirac_prod_spinspin(temp,m,in);
    spinspin_copy(out,temp);
  }
  
  //out+=m*in
  inline void spinspin_dirac_summ_the_prod_spinspin(spinspin out,const dirac_matr& in,const spinspin c)
  {
    spinspin temp;
    unsafe_dirac_prod_spinspin(temp,in,c);
    spinspin_summassign(out,temp);
  }
  
  //out=m*in^t
  inline void unsafe_dirac_prod_spinspin_transp(spinspin out,const dirac_matr& m,const spinspin in)
  {
    for(int id1=0;id1<NDIRAC;id1++)
      for(int id2=0;id2<NDIRAC;id2++)
	unsafe_complex_prod(out[id1][id2],m.entr[id1],in[id2][m.pos[id1]]);
  }
  inline void safe_dirac_prod_spinspin_transp(spinspin out,const dirac_matr& m,const spinspin in)
  {spinspin temp;unsafe_dirac_prod_spinspin_transp(temp,m,in);spinspin_copy(out,temp);}
  
  //out=m*in^+
  CUDA_HOST_AND_DEVICE inline void unsafe_dirac_prod_spinspin_dag(spinspin out,const dirac_matr& m,const spinspin in)
  {
    for(int id1=0;id1<NDIRAC;id1++)
      for(int id2=0;id2<NDIRAC;id2++)
	unsafe_complex_conj2_prod(out[id1][id2],m.entr[id1],in[id2][m.pos[id1]]);
  }
  inline void safe_dirac_prod_spinspin_dag(spinspin out,const dirac_matr& m,const spinspin in)
  {spinspin temp;unsafe_dirac_prod_spinspin_dag(temp,m,in);spinspin_copy(out,temp);}
  
  //out=in*m
  CUDA_HOST_AND_DEVICE inline void unsafe_spinspin_prod_dirac(spinspin out,const spinspin in,const dirac_matr& m)
  {
    spinspin_put_to_zero(out);
    for(int id1=0;id1<NDIRAC;id1++)
      for(int id2=0;id2<NDIRAC;id2++)
	unsafe_complex_prod(out[id1][m.pos[id2]],in[id1][id2],m.entr[id2]);
  }
  inline void safe_spinspin_prod_dirac(spinspin out,const spinspin in,const dirac_matr& m)
  {spinspin temp;unsafe_spinspin_prod_dirac(temp,in,m);spinspin_copy(out,temp);}
  
  //product of spinspin and spin
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
   void unsafe_spinspin_prod_spin(A&& a,
				  const B& b,
				  const C& c)
  {
    spin_put_to_zero(a);
    for(int id1=0;id1<NDIRAC;id1++)
      for(int id2=0;id2<NDIRAC;id2++)
	complex_summ_the_prod(a[id1],b[id1][id2],c[id2]);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_spinspin_prod_spin(A&& a,
			       const B& b,
			       const C& c)
  {
    spin tmp;
    unsafe_spinspin_prod_spin(tmp,b,c);
    spin_copy(a,tmp);
  }
  
  //product of spinspin and spin
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void unsafe_spin_prod_spinspin(A&& a,
				 const B& b,
				 const C& c)
  {
    spin_put_to_zero(a);
    for(int id1=0;id1<NDIRAC;id1++)
      for(int id2=0;id2<NDIRAC;id2++)
	complex_summ_the_prod(a[id1],b[id2],c[id2][id1]);
  }
  
  template <typename A,
	    typename B,
	    typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void safe_spin_prod_spinspin(A&& a,
			       const B& b,
			       const C& c)
  {
    spin tmp;
    unsafe_spin_prod_spinspin(tmp,b,c);
    spin_copy(a,tmp);
  }
  
  //Trace of the product of two spinspins
  CUDA_HOST_AND_DEVICE inline void summ_the_trace_prod_spinspins(complex c,const spinspin a,const spinspin b)
  {
    for(int id1=0;id1<NDIRAC;id1++)
      for(int id2=0;id2<NDIRAC;id2++)
	complex_summ_the_prod(c,a[id1][id2],b[id2][id1]);
  }
  CUDA_HOST_AND_DEVICE inline void trace_prod_spinspins(complex c,const spinspin a,const spinspin b)
  {
    c[0]=c[1]=0;
    summ_the_trace_prod_spinspins(c,a,b);
  }
  CUDA_HOST_AND_DEVICE inline void trace_spinspin_with_dirac(complex out,const spinspin s,const dirac_matr& m)
  {
    complex_put_to_zero(out);
    for(int id1=0;id1<NDIRAC;id1++)
      complex_summ_the_prod(out,m.entr[id1],s[m.pos[id1]][id1]);
  }
  
  //dirac*spin
  inline void unsafe_dirac_prod_spin(spin out,const dirac_matr& m,const spin in)
  {for(int id1=0;id1<NDIRAC;id1++) unsafe_complex_prod(out[id1],m.entr[id1],in[m.pos[id1]]);}
  inline void safe_dirac_prod_spin(spin out,const dirac_matr& m,const spin in){spin tmp;unsafe_dirac_prod_spin(tmp,m,in);spin_copy(out,tmp);}
  inline void unsafe_spin_prod_dirac(spin out,const spin in,const dirac_matr& m)
  {spin_put_to_zero(out);for(int id1=0;id1<NDIRAC;id1++) complex_summ_the_prod(out[m.pos[id1]],in[id1],m.entr[id1]);}
  inline void safe_spin_prod_spin(spin out,const spin in,const dirac_matr& m){spin tmp;unsafe_spin_prod_dirac(tmp,in,m);spin_copy(out,tmp);}
  
  //Rotate left and right by (1+-ig5)/sqrt(2)
  //We distinguish for types of rotations, ir=rsink*2+rsource
  //We divide the spinspin into 4 blocks, according to id<2
  // -------------------------------
  // | div |  0  |  1  |  2  |  3  |
  // -------------------------------
  // | 0 1 | + 1 | 1 + | 1 - | - 1 |
  // | 2 3 | 1 - | - 1 | + 1 | 1 + |
  // -------------------------------
  // so just have to specify which is the block which rotate as +-i
  CUDA_HOST_AND_DEVICE inline void rotate_spinspin_to_physical_basis(spinspin s,int rsi,int rso)
  {
    const int list_prb[4]={0,1,2,3},list_mrb[4]={3,2,1,0}; //plus and minus rotating blocks
    const int so_shft[4]={0,2,0,2},si_shft[4]={0,0,2,2};   //start of dirac indexes defining blocks
    
    int ir=rsi*2+rso,prb=list_prb[ir],mrb=list_mrb[ir];
    
    for(int dso=0;dso<2;dso++)
      for(int dsi=0;dsi<2;dsi++)
	{
	  int pso=dso+so_shft[prb],psi=dsi+si_shft[prb];
	  int mso=dso+so_shft[mrb],msi=dsi+si_shft[mrb];
	  
	  // rotation with +,-
	  assign_complex_prod_i(s[pso][psi]);
	  assign_complex_prod_minus_i(s[mso][msi]);
	}
  }
  
  //compute the determinant of spinspin matrix
  inline void spinspin_det(complex d,const spinspin s)
  {matrix_determinant(d,(complex*)s,4);}
  
  /// a=b
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spin1field_copy(A&& a,
		       const B& b)
  {
    for(int mu=0;mu<NDIM;mu++)
      complex_copy(a[mu],b[mu]);
  }
  
  /// a=b
  template <typename A,
	    typename B>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void spin1prop_copy(A&& a,
		      const B& b)
  {
    for(int mu=0;mu<NDIM;mu++)
      spin1field_copy(a[mu],b[mu]);
  }
}

#undef EXTERN_SPIN

#endif
