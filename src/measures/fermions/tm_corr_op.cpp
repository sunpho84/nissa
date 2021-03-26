#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "linalgs/reduce.hpp"

#include "tm_corr_op.hpp"

namespace nissa
{
  void tm_corr_op::ins(spincolor *out,const int igamma,spincolor *in)
  {
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	unsafe_dirac_prod_spincolor(out[ivol],base_gamma+igamma,in[ivol]);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  void tm_corr_op::undiluted_meson_contr(complex* contr,
					 spincolor *bw,
					 spincolor *fw,
					 const int& igamma,
					 const int source_coord)
  {
    /// Local storage
    complex *loc_contr=get_reducing_buffer<complex>(loc_vol);
    
    /// Gamma matrix
    dirac_matr g=base_gamma[igamma]*base_gamma[5];
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	spincolor temp;
	unsafe_dirac_prod_spincolor(temp,&g,fw[ivol]);
	spincolor_scalar_prod(loc_contr[ivol],bw[ivol],fw[ivol]);
      }
    NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
    
    complex unshifted_glb_contr[glb_size[0]];
    glb_reduce(unshifted_glb_contr,loc_contr,loc_vol,glb_size[0],loc_size[0],glb_coord_of_loclx[0][0]);
    
    for(int glb_t=0;glb_t<glb_size[0];glb_t++)
      {
	/// Distance from source
	const int dt=
	  (glb_t-source_coord+glb_size[0])%glb_size[0];
	
	complex_copy(contr[dt],unshifted_glb_contr[glb_t]);
      }
    
  }
}

