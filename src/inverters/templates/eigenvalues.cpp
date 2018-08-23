#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "eigenvalues.hpp"

//part of this should be moved to linalgs?

namespace nissa
{
  namespace internal_eigenvalues
  {
    //scalar product of (in1,in2)
    void scalar_prod(complex out,complex *in1,complex *in2,complex *buffer,int loc_size)
    {
      GET_THREAD_ID();
      
      //local scalar product
      NISSA_PARALLEL_LOOP(i,0,loc_size)
	unsafe_complex_conj1_prod(buffer[i],in1[i],in2[i]);
      THREAD_BARRIER();
      
      //reduction
      complex_vector_glb_collapse(out,buffer,loc_size);
    }
    
    //Ai=Ai-Bi*c
    void complex_vector_subtssign_complex_vector_prod_complex(complex *a,complex *b,complex c,int n)
    {
      GET_THREAD_ID();
      
      NISSA_PARALLEL_LOOP(i,0,n)
	complex_subt_the_prod(a[i],b[i],c);
      set_borders_invalid(a);
    }
    
    //orthogonalize v with respect to the nvec of V
    void modified_GS(complex *v,complex **V,complex *buffer,int nvec,int vec_size)
    {
      for(int i=0;i<nvec;i++)
	{
	  complex s;
	  scalar_prod(s,V[i],v,buffer,vec_size);
	  complex_vector_subtssign_complex_vector_prod_complex(v,V[i],s,vec_size);
	}
    }
  }
}
